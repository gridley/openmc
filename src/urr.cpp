#include "openmc/urr.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/simulation.h"

#include <algorithm> // any_of
#include <iostream>

namespace openmc {

UrrData::UrrData(hid_t group_id)
{
  // Read interpolation and other flags
  int interp_temp;
  read_attribute(group_id, "interpolation", interp_temp);
  interp_ = static_cast<Interpolation>(interp_temp);

  // read the metadata
  read_attribute(group_id, "inelastic", inelastic_flag_);
  read_attribute(group_id, "absorption", absorption_flag_);
  int temp_multiply_smooth;
  read_attribute(group_id, "multiply_smooth", temp_multiply_smooth);
  multiply_smooth_ = (temp_multiply_smooth == 1);

  // read the energies at which tables exist
  read_dataset(group_id, "energy", energy_);

  // Read URR tables. The HDF5 format is a little
  // different from how we want it laid out in memory.
  // This array used to be called "prob_".
  xt::xtensor<double, 3> tmp_prob;
  read_dataset(group_id, "table", tmp_prob);
  auto shape = tmp_prob.shape();

  // We separate out into two matrices (one with CDF values,
  // the other with cross section sets) in order to improve
  // contiguity of memory accesses.
  const auto n_energy = shape[0];
  const auto n_cdf_values = shape[2];
  cdf_values_.resize({n_energy, n_cdf_values});
  xs_values_.resize({n_energy, n_cdf_values});

  // Now fill in the values. Using manual loops here since we might
  // not have fancy xtensor slicing code written for GPU tensors.
  // The below enum gives how URR tables are laid out in our HDF5 tables.
  enum class URRTableParam {
    CUM_PROB,
    TOTAL,
    ELASTIC,
    FISSION,
    N_GAMMA,
    HEATING
  };
  for (int energy_index = 0; energy_index < n_energy; ++energy_index) {
    for (int i_cdf = 0; i_cdf < n_cdf_values; ++i_cdf) {
      cdf_values_(energy_index, i_cdf) =
        tmp_prob(energy_index, URRTableParam::CUM_PROB, i_cdf);
      xs_values_(energy_index, i_cdf).total =
        tmp_prob(energy_index, URRTableParam::TOTAL, i_cdf);
      xs_values_(energy_index, i_cdf).elastic =
        tmp_prob(energy_index, URRTableParam::ELASTIC, i_cdf);
      xs_values_(energy_index, i_cdf).fission =
        tmp_prob(energy_index, URRTableParam::FISSION, i_cdf);
      xs_values_(energy_index, i_cdf).n_gamma =
        tmp_prob(energy_index, URRTableParam::N_GAMMA, i_cdf);
      xs_values_(energy_index, i_cdf).heating =
        tmp_prob(energy_index, URRTableParam::HEATING, i_cdf);
    }
  }
}

bool UrrData::has_negative() const
{

  // Lambda checks if any value in XSSset is negative
  auto xs_set_negative = [](const XSSet& xs) {
    return xs.total < 0.0 || xs.elastic < 0.0 || xs.fission < 0.0 ||
           xs.n_gamma < 0.0 || xs.heating < 0.0;
  };

  return std::any_of(cdf_values_.begin(), cdf_values_.end(), [](double x) {
    return x < 0.0;
  }) || std::any_of(xs_values_.begin(), xs_values_.end(), xs_set_negative);
}

// Samples an inverse Gaussian random variable
double sample_ig(double mu, double lam, uint64_t* seed) {

  // This implementation of normal_variate takes an
  // undefined number of prn() calls, so the URR state
  // is fast-forwarded assuming around 100 times. There
  // might be VERY small correlations, but not likely anything
  // that matters for particle transport.
  double w = mu * std::pow(normal_variate(0.0, 1.0, seed), 2);
  double c = 0.5 * mu / lam;
  double x1 = mu + c * (w - std::sqrt(w*(4*lam+w)));
  double x = x1;
  if (prn(seed) >= mu / (mu + x1)) {
    x = mu * mu / x1;
  }
  return x;
}

// Samples a normal inverse Gaussian random variable
double sample_nig(double alpha, double beta, double mu, double delta2, uint64_t* seed) {
  assert(std::abs(alpha) > std::abs(beta));
  double z = sample_ig(std::sqrt(delta2 / (alpha*alpha-beta*beta)), delta2, seed);
  return std::sqrt(z) * normal_variate(0.0, 1.0, seed) + beta * z + mu;
}

ContinuousURRData::ContinuousURRData(const std::string& filename, gsl::index index) {
  hid_t h5file = file_open(filename.c_str(), 'r', false);
  hid_t h5curvefit = open_group(h5file, "total_xs_parameters");
  hid_t h5conditionals = open_group(h5file, "conditional_partials");
  read_dataset(h5curvefit, "alpha", alpha, true);
  read_dataset(h5curvefit, "beta", beta, true);
  read_dataset(h5curvefit, "mu", mu, true);
  read_dataset(h5curvefit, "delta2", delta2, true);
  read_dataset(h5conditionals, "x_values", nodes, true);
  read_dataset(h5conditionals, "w_values", weights, true);
  read_dataset(h5conditionals, "absorption", abs_values, true);
  if (object_exists(h5conditionals, "fission")) {
    read_dataset(h5conditionals, "fission", fiss_values, true);
    has_fission_ = true;
  }

  read_dataset(h5file, "energies", energy_, true);

  file_close(h5file);


  // Save the nuclide index to keep the same LCG stream between lookups
  // at this energy.
  index_ = index;
}

void ContinuousURRData::sample(double E, int i_T, uint64_t* seed, NuclideMicroXS& xs) {
  // Look up the energy index to use here
  if (!energy_in_bounds(E)) {
    fatal_error("energy out of bounds in continuous URR");
  }

  int energy_index = lower_bound_index(energy_.begin(), energy_.end(), E);

  // We use more random numbers than in the single table case. It's
  // an undefined number but quite likely less than 400.
  uint64_t fseed = future_seed(static_cast<uint64_t>(400 * index_), *seed);

  // energy interpolation factor
  double f = (E - energy_[energy_index]) /
        (energy_[energy_index + 1] - energy_[energy_index]);

  double a = (1.0 - f) * alpha(energy_index, i_T) + f * alpha(energy_index + 1, i_T);
  double b = (1.0 - f) * beta(energy_index, i_T) + f * beta(energy_index + 1, i_T);
  double m = (1.0 - f) * mu(energy_index, i_T) + f * mu(energy_index + 1, i_T);
  double d2 =(1.0 - f) * delta2(energy_index, i_T) + f * delta2(energy_index + 1, i_T);

  double sigt = sample_nig(0.5 * (a + b), -0.5 * (b - a), m, d2, &fseed);

  if (sigt <= 0.0) {
    sigt = 1e-4;
    xs.total = 1e-4;
    xs.absorption = 1e-4;
    xs.fission = 0.0;
    xs.nu_fission = 0.0;
    xs.elastic = 0.0;
    return;
  }

  // Note that "abs" has been used to represent "capture".

  constexpr int bary_order = 5;
  double denom = 0.0;
  double abs = 0.0;
  double fiss = 0.0;

  for (int j=0; j<bary_order; ++j) {
    double term = weights(energy_index, j) / (sigt - nodes(energy_index, j));
    abs += abs_values(energy_index, j, i_T) * term;
    if (has_fission_)
      fiss += fiss_values(energy_index, j, i_T) * term;
    denom += term;
  }
  abs /= denom;
  fiss /= denom;

  double mul = uniform_distribution(0.99, 1.01, &fseed);
  xs.total *= mul;
  xs.absorption *= mul;
  xs.fission *= mul;
  xs.elastic *= mul;

  // xs.total = sigt;
  // xs.absorption = abs + fiss;
  // xs.fission = fiss;
  // xs.elastic = sigt - abs - fiss;

  // if (simulation::need_depletion_rx) {
  //   // Separate the pure capture component
  //   xs.reaction[0] = abs;
  // }
}

} // namespace openmc
