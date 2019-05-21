#include "openmc/thermal.h"

#include <algorithm> // for sort, move, min, max, find
#include <cmath>     // for round, sqrt, abs
#include <sstream>   // for stringstream

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_thermal.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {
std::vector<std::unique_ptr<ThermalScattering>> thermal_scatt;
std::unordered_map<std::string, int> thermal_scatt_map;
}

//==============================================================================
// ThermalScattering implementation
//==============================================================================

ThermalScattering::ThermalScattering(hid_t group, const std::vector<double>& temperature)
{
  // Get name of table from group
  name_ = object_name(group);

  // Get rid of leading '/'
  name_ = name_.substr(1);

  read_attribute(group, "atomic_weight_ratio", awr_);
  read_attribute(group, "nuclides", nuclides_);

  // Read temperatures
  hid_t kT_group = open_group(group, "kTs");

  // Determine temperatures available
  auto dset_names = dataset_names(kT_group);
  auto n = dset_names.size();
  auto temps_available = xt::empty<double>({n});
  for (int i = 0; i < dset_names.size(); ++i) {
    // Read temperature value
    double T;
    read_dataset(kT_group, dset_names[i].data(), T);
    temps_available[i] = T / K_BOLTZMANN;
  }
  std::sort(temps_available.begin(), temps_available.end());

  // Determine actual temperatures to read -- start by checking whether a
  // temperature range was given, in which case all temperatures in the range
  // are loaded irrespective of what temperatures actually appear in the model
  std::vector<int> temps_to_read;
  if (settings::temperature_range[1] > 0.0) {
    for (const auto& T : temps_available) {
      if (settings::temperature_range[0] <= T &&
          T <= settings::temperature_range[1]) {
        temps_to_read.push_back(std::round(T));
      }
    }
  }

  switch (settings::temperature_method) {
  case TEMPERATURE_NEAREST:
    // Determine actual temperatures to read
    for (const auto& T : temperature) {

      auto i_closest = xt::argmin(xt::abs(temps_available - T))[0];
      auto temp_actual = temps_available[i_closest];
      if (std::abs(temp_actual - T) < settings::temperature_tolerance) {
        if (std::find(temps_to_read.begin(), temps_to_read.end(), std::round(temp_actual))
            == temps_to_read.end()) {
          temps_to_read.push_back(std::round(temp_actual));
        }
      } else {
        std::stringstream msg;
        msg << "Nuclear data library does not contain cross sections for "
          << name_ << " at or near " << std::round(T) << " K.";
        fatal_error(msg);
      }
    }
    break;

  case TEMPERATURE_INTERPOLATION:
    // If temperature interpolation or multipole is selected, get a list of
    // bounding temperatures for each actual temperature present in the model
    for (const auto& T : temperature) {
      bool found = false;
      for (int j = 0; j < temps_available.size() - 1; ++j) {
        if (temps_available[j] <= T && T < temps_available[j + 1]) {
          int T_j = std::round(temps_available[j]);
          int T_j1 = std::round(temps_available[j + 1]);
          if (std::find(temps_to_read.begin(), temps_to_read.end(), T_j) == temps_to_read.end()) {
            temps_to_read.push_back(T_j);
          }
          if (std::find(temps_to_read.begin(), temps_to_read.end(), T_j1) == temps_to_read.end()) {
            temps_to_read.push_back(T_j1);
          }
          found = true;
        }
      }
      if (!found) {
        std::stringstream msg;
        msg << "Nuclear data library does not contain cross sections for "
          << name_ << " at temperatures that bound " << std::round(T) << " K.";
        fatal_error(msg);
      }
    }
  }

  // Sort temperatures to read
  std::sort(temps_to_read.begin(), temps_to_read.end());

  auto n_temperature = temps_to_read.size();
  kTs_.reserve(n_temperature);
  data_.reserve(n_temperature);

  for (auto T : temps_to_read) {
    // Get temperature as a string
    std::string temp_str = std::to_string(T) + "K";

    // Read exact temperature value
    double kT;
    read_dataset(kT_group, temp_str.data(), kT);
    kTs_.push_back(kT);

    // Open group for temperature i
    hid_t T_group = open_group(group, temp_str.data());
    data_.emplace_back(T_group);
    close_group(T_group);
  }

  close_group(kT_group);
}

void
ThermalScattering::calculate_xs(double E, double sqrtkT, int* i_temp,
                                double* elastic, double* inelastic) const
{
  // Determine temperature for S(a,b) table
  double kT = sqrtkT*sqrtkT;
  int i;
  if (settings::temperature_method == TEMPERATURE_NEAREST) {
    // If using nearest temperature, do linear search on temperature
    for (i = 0; i < kTs_.size(); ++i) {
      if (std::abs(kTs_[i] - kT) < K_BOLTZMANN*settings::temperature_tolerance) {
        break;
      }
    }
  } else {
    // Find temperatures that bound the actual temperature
    for (i = 0; i < kTs_.size() - 1; ++i) {
      if (kTs_[i] <= kT && kT < kTs_[i+1]) {
        break;
      }
    }

    // Randomly sample between temperature i and i+1
    double f = (kT - kTs_[i]) / (kTs_[i+1] - kTs_[i]);
    if (f > prn()) ++i;
  }

  // Set temperature index
  *i_temp = i;

  // Get pointer to S(a,b) table
  auto& sab = data_[i];

  // Get index and interpolation factor for inelastic grid
  int i_grid = 0;
  double f = 0.0;
  if (E >= sab.inelastic_e_in_.front()) {
    auto& E_in = sab.inelastic_e_in_;
    i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
    f = (E - E_in[i_grid]) / (E_in[i_grid+1] - E_in[i_grid]);
  }

  // Calculate S(a,b) inelastic scattering cross section
  auto& xs = sab.inelastic_sigma_;
  *inelastic = xs[i_grid] + f * (xs[i_grid + 1] - xs[i_grid]);

  // Check for elastic data
  if (!sab.elastic_e_in_.empty()) {
    // Determine whether elastic scattering is given in the coherent or
    // incoherent approximation. For coherent, the cross section is
    // represented as P/E whereas for incoherent, it is simply P

    auto& E_in = sab.elastic_e_in_;

    if (sab.elastic_mode_ == SAB_ELASTIC_COHERENT) {
      if (E < E_in.front()) {
        // If energy is below that of the lowest Bragg peak, the elastic
        // cross section will be zero
        *elastic = 0.0;
      } else {
        i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
        *elastic = sab.elastic_P_[i_grid] / E;
      }
    } else {
      // Determine index on elastic energy grid
      if (E < E_in.front()) {
        i_grid = 0;
      } else {
        i_grid = lower_bound_index(E_in.begin(), E_in.end(), E);
      }

      // Get interpolation factor for elastic grid
      f = (E - E_in[i_grid]) / (E_in[i_grid+1] - E_in[i_grid]);

      // Calculate S(a,b) elastic scattering cross section
      auto& xs = sab.elastic_P_;
      *elastic = xs[i_grid] + f*(xs[i_grid + 1] - xs[i_grid]);
    }
  } else {
    // No elastic data
    *elastic = 0.0;
  }
}

bool
ThermalScattering::has_nuclide(const char* name) const
{
  std::string nuc {name};
  return std::find(nuclides_.begin(), nuclides_.end(), nuc) != nuclides_.end();
}

//==============================================================================
// ThermalData implementation
//==============================================================================

ThermalData::ThermalData(hid_t group)
{
  // Coherent/incoherent elastic data
  if (object_exists(group, "elastic")) {
    // Read cross section data
    hid_t elastic_group = open_group(group, "elastic");

    // Read elastic cross section
    elastic_.xs = read_function(elastic_group, "xs");

    // Read angle-energy distribution
    hid_t dgroup = open_group(elastic_group, "distribution");
    std::string temp;
    read_attribute(dgroup, "type", temp);
    if (temp == "coherent_elastic") {
      auto xs = dynamic_cast<CoherentElasticXS*>(elastic_.xs.get());
      elastic_.distribution = std::make_unique<CoherentElasticAE>(*xs);

      // Set threshold energy
      threshold_elastic_ = xs->bragg_edges().back();

    } else {
      auto xs = dynamic_cast<Tabulated1D*>(elastic_.xs.get());
      if (temp == "incoherent_elastic") {
        elastic_.distribution = std::make_unique<IncoherentElasticAE>(dgroup);
      } else if (temp == "incoherent_elastic_discrete") {
        elastic_.distribution = std::make_unique<IncoherentElasticAEDiscrete>(
          dgroup, xs->x()
        );
      }

      // Set threshold energy
      threshold_elastic_ = xs->x().back();
    }

    close_group(elastic_group);
  }

  // Inelastic data
  if (object_exists(group, "inelastic")) {
    // Read type of inelastic data
    hid_t inelastic_group = open_group(group, "inelastic");

    // Read inelastic cross section
    inelastic_.xs = read_function(inelastic_group, "xs");

    // Set inelastic threshold
    auto xs = dynamic_cast<Tabulated1D*>(inelastic_.xs.get());
    threshold_inelastic_ = xs->x().back();

    // Read angle-energy distribution
    hid_t dgroup = open_group(inelastic_group, "distribution");
    std::string temp;
    read_attribute(dgroup, "type", temp);
    if (temp == "incoherent_inelastic") {
      inelastic_.distribution = std::make_unique<IncoherentInelasticAE>(dgroup);
    } else if (temp == "incoherent_inelastic_discrete") {
      inelastic_.distribution = std::make_unique<IncoherentInelasticAEDiscrete>(
        dgroup, xs->x()
      );
    }

    close_group(inelastic_group);
  }
}

void
ThermalData::sample(const NuclideMicroXS& micro_xs, double E,
                    double* E_out, double* mu)
{
  // Determine whether inelastic or elastic scattering will occur
  if (prn() < micro_xs.thermal_elastic / micro_xs.thermal) {
    elastic_.distribution->sample(E, *E_out, *mu);
  } else {
    inelastic_.distribution->sample(E, *E_out, *mu);
  }

  // Because of floating-point roundoff, it may be possible for mu to be
  // outside of the range [-1,1). In these cases, we just set mu to exactly
  // -1 or 1
  if (std::abs(*mu) > 1.0) *mu = std::copysign(1.0, *mu);
}

void free_memory_thermal()
{
  data::thermal_scatt.clear();
  data::thermal_scatt_map.clear();
}

} // namespace openmc
