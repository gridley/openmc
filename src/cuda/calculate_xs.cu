
#include <thrust/complex.h>

#include "openmc/cuda/calculate_xs.h"
#include "openmc/geometry.h"         // find_cell
#include "openmc/reaction_product.h" // EmissionMode
#include "openmc/search.h"
#include "openmc/settings.h" // BLOCKSIZE

#include "openmc/wmp.h" // TODO remove?

// TODO potentially load p.sqrtkT() to a register

__device__ thrust::complex<double> zpf8h_faddeeva(thrust::complex<double> z)
{
  double flip_real_part =
    1.0; // TODO test sign flip with bit magic. Just set mask based on..
  if (z.imag() < 0.0) {
    flip_real_part = -1.0;
    z.imag(-z.imag()); // bit magic?
  }
  z.imag(z.imag() + 0.9);
  const auto zz = z * z;

  const double z_r = z.real();
  const double z_i = z.imag();
  const double zz_r = zz.real();
  const double zz_i = zz.imag();

  constexpr double aa0_r = 11.7559071436993;
  constexpr double aa1_i = -32.310199761603;
  constexpr double aa2_r = -21.9357456686406;
  constexpr double aa3_i = 31.490536152863;
  constexpr double aa4_r = 6.75847413957232;
  constexpr double aa5_i = -8.07354660639634;
  constexpr double aa6_r = -0.507771291744591;
  constexpr double aa7_i = 0.564189504758109;

  constexpr double bb0_r = 6.5625;
  constexpr double bb1_r = -52.5;
  constexpr double bb2_r = 52.5;
  constexpr double bb3_r = -14.0;

  const double num_re =
    (((((((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
             ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
               z_i) +
            aa4_r) *
             z_r -
           ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
             ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
               z_r) *
             z_i) *
            z_r -
          (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_i) +
              aa4_r) *
               z_i +
             ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_r) *
               z_r) +
            aa3_i) *
            z_i) +
         aa2_r) *
          z_r -
        (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
             ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
               z_i) +
            aa4_r) *
             z_r -
           ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
             ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
               z_r) *
             z_i) *
            z_i +
          (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_i) +
              aa4_r) *
               z_i +
             ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_r) *
               z_r) +
            aa3_i) *
            z_r) *
          z_i) *
         z_r -
       ((((((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_i) +
              aa4_r) *
               z_r -
             ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_r) *
               z_i) *
              z_r -
            (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
                 ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                   aa5_i) *
                   z_i) +
                aa4_r) *
                 z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
                 ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                   aa5_i) *
                   z_r) *
                 z_r) +
              aa3_i) *
              z_i) +
           aa2_r) *
            z_i +
          (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_i) +
              aa4_r) *
               z_r -
             ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                 aa5_i) *
                 z_r) *
               z_i) *
              z_i +
            (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
                 ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                   aa5_i) *
                   z_i) +
                aa4_r) *
                 z_i +
               ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
                 ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                   aa5_i) *
                   z_r) *
                 z_r) +
              aa3_i) *
              z_r) *
            z_r) +
         aa1_i) *
         z_i) +
      aa0_r);
  const double num_im =
    ((((((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
            ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
              z_i) +
           aa4_r) *
            z_r -
          ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
            ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
              z_r) *
            z_i) *
           z_r -
         (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_i) +
             aa4_r) *
              z_i +
            ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_r) *
              z_r) +
           aa3_i) *
           z_i) +
        aa2_r) *
         z_r -
       (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
            ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
              z_i) +
           aa4_r) *
            z_r -
          ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
            ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
              z_r) *
            z_i) *
           z_i +
         (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_i) +
             aa4_r) *
              z_i +
            ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_r) *
              z_r) +
           aa3_i) *
           z_r) *
         z_i) *
        z_i +
      ((((((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_i) +
             aa4_r) *
              z_r -
            ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_r) *
              z_i) *
             z_r -
           (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
                ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                  aa5_i) *
                  z_i) +
               aa4_r) *
                z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
                ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                  aa5_i) *
                  z_r) *
                z_r) +
             aa3_i) *
             z_i) +
          aa2_r) *
           z_i +
         (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_i) +
             aa4_r) *
              z_r -
            ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) + aa5_i) *
                z_r) *
              z_i) *
             z_i +
           (((((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_r -
                ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                  aa5_i) *
                  z_i) +
               aa4_r) *
                z_i +
              ((((-aa7_i * z_i) + aa6_r) * z_r - (aa7_i * z_r) * z_i) * z_i +
                ((((-aa7_i * z_i) + aa6_r) * z_i + (aa7_i * z_r) * z_r) +
                  aa5_i) *
                  z_r) *
                z_r) +
             aa3_i) *
             z_r) *
           z_r) +
        aa1_i) *
        z_r);
  const double den_re =
    (((((((zz_r + bb3_r) * zz_r - zz_i * zz_i) + bb2_r) * zz_r -
         ((zz_r + bb3_r) * zz_i + zz_i * zz_r) * zz_i) +
        bb1_r) *
         zz_r -
       ((((zz_r + bb3_r) * zz_r - zz_i * zz_i) + bb2_r) * zz_i +
         ((zz_r + bb3_r) * zz_i + zz_i * zz_r) * zz_r) *
         zz_i) +
      bb0_r);
  const double den_im =
    ((((((zz_r + bb3_r) * zz_r - zz_i * zz_i) + bb2_r) * zz_r -
        ((zz_r + bb3_r) * zz_i + zz_i * zz_r) * zz_i) +
       bb1_r) *
        zz_i +
      ((((zz_r + bb3_r) * zz_r - zz_i * zz_i) + bb2_r) * zz_i +
        ((zz_r + bb3_r) * zz_i + zz_i * zz_r) * zz_r) *
        zz_r);
  const double modulus = den_re * den_re + den_im * den_im;
  return {flip_real_part * (num_re * den_re + num_im * den_im) / modulus,
    (num_im * den_re - num_re * den_im) / modulus};
}

namespace openmc {
namespace gpu {

__constant__ unique_ptr<Material>* materials;
__constant__ unique_ptr<Nuclide>* nuclides;
__constant__ unique_ptr<ThermalScattering>* thermal_scatt;
__constant__ Particle* particles;
__constant__ xsfloat energy_min_neutron;
__constant__ xsfloat energy_max_neutron;
__constant__ xsfloat log_spacing;
__constant__ unsigned number_nuclides;
__constant__ bool need_depletion_rx;

__managed__ unsigned managed_calculate_fuel_queue_index;
__managed__ unsigned managed_calculate_nonfuel_queue_index;

__global__ void __launch_bounds__(BLOCKSIZE) process_calculate_xs_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size)
{
  using EmissionMode = ReactionProduct::EmissionMode;

  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  if (tid >= queue_size)
    return;
  Particle p(queue[tid].idx);
  auto const E = __ldg(&queue[tid].E);
  auto const mat_idx = __ldg(&queue[tid].material);

  // Store pre-collision particle properties
  p.wgt_last() = p.wgt();
  p.E_last() = E;
  p.u_last() = p.u();
  p.r_last() = p.r();

  // Reset event variables
  p.event() = TallyEvent::KILL;
  p.event_nuclide() = NUCLIDE_NONE;
  p.event_mt() = REACTION_NONE;

  p.macro_xs().total = 0.0;
  p.macro_xs().neutron.absorption = 0.0;
  p.macro_xs().neutron.fission = 0.0;
  p.macro_xs().neutron.nu_fission = 0.0;

  // Skip void material
  if (mat_idx == -1)
    return;

  Material const& m = *materials[mat_idx];

  // Used for looping over thermal tables. Tryna use less registers
  unsigned char n_therm_tables = m.thermal_tables_.size();
  unsigned char therm_table_indx = 0;
  int next_sab_nuclide = -1;
  if (n_therm_tables) {
    next_sab_nuclide = m.thermal_tables_[therm_table_indx].index_nuclide;
  }

  unsigned i_log_union = std::log(E / energy_min_neutron) / log_spacing;

  // Add contribution from each nuclide in material
  auto const n_nuclides = m.nuclide_.size();
  for (int i = 0; i < n_nuclides; ++i) {

    auto const& i_nuclide =
      m.nuclide_[i]; // TODO test if making not a reference better
    auto* __restrict__ micro_ref {&p.neutron_xs(i_nuclide)};
    NuclideMicroXS micro;

    micro.index_sab = C_NONE;
    micro.sab_frac = 0.0;
    if (i == next_sab_nuclide) {
      const auto sab {m.thermal_tables_[therm_table_indx]};
      micro.index_sab = sab.index_table;
      micro.sab_frac = sab.fraction;
      if (E > gpu::thermal_scatt[micro.index_sab]->energy_max_)
        micro.index_sab = C_NONE;
      ++therm_table_indx;
      if (therm_table_indx < n_therm_tables) {
        next_sab_nuclide = m.thermal_tables_[therm_table_indx].index_nuclide;
      } else {
        next_sab_nuclide = -1; // done with S(a, b)
      }
    }

    NuclideMicroXS* __restrict__ use_micro = micro_ref;
    if (E != micro_ref->last_E || p.sqrtkT() != micro_ref->last_sqrtkT ||
        micro.index_sab != micro_ref->index_sab ||
        micro.sab_frac != micro_ref->sab_frac) {
      use_micro = &micro; // ummm
      auto const& nuclide = *nuclides[i_nuclide];
      micro.elastic = CACHE_INVALID;
      micro.thermal = 0.0;
      micro.thermal_elastic = 0.0;
      micro.use_ptable = false;
      micro.last_E = E;
      micro.last_sqrtkT = p.sqrtkT();

      // Find the appropriate temperature index. why would someone use
      // nearest?
      xsfloat kT = p.sqrtkT() * p.sqrtkT();
      xsfloat f;
      int i_temp;

      switch (gpu::temperature_method) {
      case TemperatureMethod::NEAREST: {
        double max_diff = INFTY;
        for (int t = 0; t < nuclide.kTs_.size(); ++t) {
          double diff = std::abs(nuclide.kTs_[t] - kT);
          if (diff < max_diff) {
            i_temp = t;
            max_diff = diff;
          }
        }
      } break;

      case TemperatureMethod::INTERPOLATION:
        // Find temperatures that bound the actual temperature
        for (i_temp = 0; i_temp < nuclide.kTs_.size() - 1; ++i_temp) {
          if (nuclide.kTs_[i_temp] <= kT && kT < nuclide.kTs_[i_temp + 1])
            break;
        }

        // Randomly sample between temperature i and i+1
        f = (kT - nuclide.kTs_[i_temp]) /
            (nuclide.kTs_[i_temp + 1] - nuclide.kTs_[i_temp]);
        if (f > prn(p.current_seed()))
          ++i_temp;
        break;
      }

      const auto& grid {nuclide.grid_[i_temp]};
      // Determine bounding indices based on which equal log-spaced
      // interval the energy is in
      int i_low = __ldg(&grid.grid_index[i_log_union]);
      int i_high = __ldg(&grid.grid_index[i_log_union + 1]) + 1;

      // Perform binary search over reduced range
      int i_grid = i_low + lower_bound_index_linear(
                             &grid.energy[i_low], &grid.energy[i_high], E);
      const auto xs_left {nuclide.xs_[i_temp][i_grid]};
      const auto xs_right {nuclide.xs_[i_temp][i_grid + 1]};
      // check for rare case where two energy points are the same
      if (grid.energy[i_grid] == grid.energy[i_grid + 1])
        ++i_grid;

      // calculate interpolation factor
      f = (E - grid.energy[i_grid]) /
          (grid.energy[i_grid + 1] - grid.energy[i_grid]);

      micro.index_temp = i_temp;
      micro.index_grid = i_grid;
      micro.interp_factor = f;

      // Calculate all microscopic cross sections
      micro.total = (1.0 - f) * xs_left.total + f * xs_right.total;
      micro.absorption =
        (1.0 - f) * xs_left.absorption + f * xs_right.absorption;

      if (nuclide.fissionable_) {
        // Calculate microscopic nuclide total cross section
        micro.fission = (1.0 - f) * xs_left.fission + f * xs_right.fission;

        // Calculate microscopic nuclide nu-fission cross section
        micro.nu_fission =
          (1.0 - f) * xs_left.nu_fission + f * xs_right.nu_fission;
      } else {
        micro.fission = 0.0;
        micro.nu_fission = 0.0;
      }

      // Calculate microscopic nuclide photon production cross section
      micro.photon_prod =
        (1.0 - f) * xs_left.photon_production + f * xs_right.photon_production;

      // Additionally calculate S(a, b) cross section data
      if (micro.index_sab >= 0) {
        int i_temp;
        xsfloat inelastic;
        xsfloat elastic;
        gpu::thermal_scatt[micro.index_sab]->calculate_xs(E, micro.last_sqrtkT,
          &i_temp, &elastic, &inelastic, p.current_seed());
        micro.thermal = micro.sab_frac * (elastic + inelastic);
        micro.thermal_elastic = micro.sab_frac * elastic;

        // calculate_elastic_xs
        if (micro.index_temp >= 0) {
          const auto& xs = nuclide.reactions_[0]->xs_[micro.index_temp].value;
          micro.elastic = (1.0 - micro.interp_factor) * xs[micro.index_grid] +
                          micro.interp_factor * xs[micro.index_grid + 1];
        }

        micro.total =
          micro.total + micro.thermal - micro.sab_frac * micro.elastic;
        micro.elastic = micro.thermal + (1.0 - micro.sab_frac) * micro.elastic;
        micro.index_temp_sab = i_temp;
      }

      // Calculate URR cross sections if needed
      if (gpu::urr_ptables_on && nuclide.urr_present_) {
        if (nuclide.urr_data_[i_temp].energy_in_bounds(E)) {
          micro.use_ptable = true;
          // TODO check storing by value
          const auto& urr = nuclide.urr_data_[i_temp];

          int i_energy = 0;
          while (E >= urr.energy_[i_energy + 1]) {
            ++i_energy;
          };

          p.stream() = STREAM_URR_PTABLE;
          double r =
            future_prn(static_cast<int64_t>(nuclide.index_), *p.current_seed());
          p.stream() = STREAM_TRACKING;

          int i_low = 0;
          while (urr.cdf_values_(i_energy, i_low) <= r) {
            ++i_low;
          };

          int i_up = 0;
          while (urr.cdf_values_(i_energy + 1, i_up) <= r) {
            ++i_up;
          };

          // Determine elastic, fission, and capture cross sections from the
          // probability table
          xsfloat elastic = 0.;
          xsfloat fission = 0.;
          xsfloat capture = 0.;
          xsfloat f;
          if (urr.interp_ == Interpolation::lin_lin) {
            // Determine the interpolation factor on the table
            f = (E - urr.energy_[i_energy]) /
                (urr.energy_[i_energy + 1] - urr.energy_[i_energy]);

            elastic = (1. - f) * urr.xs_values_(i_energy, i_low).elastic +
                      f * urr.xs_values_(i_energy + 1, i_up).elastic;
            fission = (1. - f) * urr.xs_values_(i_energy, i_low).fission +
                      f * urr.xs_values_(i_energy + 1, i_up).fission;
            capture = (1. - f) * urr.xs_values_(i_energy, i_low).n_gamma +
                      f * urr.xs_values_(i_energy + 1, i_up).n_gamma;
          } else if (urr.interp_ == Interpolation::log_log) {
            // Determine interpolation factor on the table
            f = std::log(E / urr.energy_[i_energy]) /
                std::log(urr.energy_[i_energy + 1] / urr.energy_[i_energy]);

            // Calculate the elastic cross section/factor
            if ((urr.xs_values_(i_energy, i_low).elastic > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).elastic > 0.)) {
              elastic = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).elastic) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).elastic));
            } else {
              elastic = 0.;
            }

            // Calculate the fission cross section/factor
            if ((urr.xs_values_(i_energy, i_low).fission > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).fission > 0.)) {
              fission = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).fission) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).fission));
            } else {
              fission = 0.;
            }

            // Calculate the capture cross section/factor
            if ((urr.xs_values_(i_energy, i_low).n_gamma > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).n_gamma > 0.)) {
              capture = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).n_gamma) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).n_gamma));
            } else {
              capture = 0.;
            }
          }

          // Determine the treatment of inelastic scattering
          xsfloat inelastic = 0.;
          if (urr.inelastic_flag_ != C_NONE) {
            // get interpolation factor
            f = micro.interp_factor;

            // Determine inelastic scattering cross section
            Reaction* rx = nuclide.reactions_[nuclide.urr_inelastic_].get();
            int xs_index = micro.index_grid - rx->xs_[i_temp].threshold;
            if (xs_index >= 0) {
              inelastic = (1. - f) * rx->xs_[i_temp].value[xs_index] +
                          f * rx->xs_[i_temp].value[xs_index + 1];
            }
          }

          // Multiply by smooth cross-section if needed
          if (urr.multiply_smooth_) {
            const auto& xs = nuclide.reactions_[0]->xs_[i_temp].value;
            f = micro.interp_factor;
            micro.elastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];
            elastic *= micro.elastic;
            capture *= (micro.absorption - micro.fission);
            fission *= micro.fission;
          }

          // Check for negative values
          if (elastic < 0.) {
            elastic = 0.;
          }
          if (fission < 0.) {
            fission = 0.;
          }
          if (capture < 0.) {
            capture = 0.;
          }

          // Set elastic, absorption, fission, total, and capture x/s. Note that
          // the total x/s is calculated as a sum of partials instead of the
          // table-provided value
          micro.elastic = elastic;
          micro.absorption = capture + fission;
          micro.fission = fission;
          micro.total = elastic + inelastic + capture + fission;

          // if (simulation::need_depletion_rx) {
          //   micro.reaction[0] = capture;
          // }

          // Determine nu-fission cross-section
          if (nuclide.fissionable_) {
            micro.nu_fission =
              nuclide.nu(E, EmissionMode::total) * micro.fission;
          }
        }
      }
    }

    double const& atom_density = m.atom_density_[i];
    p.macro_xs().total += atom_density * use_micro->total;
    p.macro_xs().neutron.absorption += atom_density * use_micro->absorption;
    p.macro_xs().neutron.fission += atom_density * use_micro->fission;
    p.macro_xs().neutron.nu_fission += atom_density * use_micro->nu_fission;
    if (use_micro != micro_ref) {
      *micro_ref = micro; // save stack variable back to global memory
    }
  }
}

__global__ void __launch_bounds__(BLOCKSIZE) process_calculate_xs_events_device_wmp(
  EventQueueItem* __restrict__ queue, unsigned queue_size)
{
  using EmissionMode = ReactionProduct::EmissionMode;

  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  if (tid >= queue_size)
    return;
  Particle p(queue[tid].idx);
  auto const E = __ldg(&queue[tid].E);
  auto const mat_idx = __ldg(&queue[tid].material);

  // Store pre-collision particle properties
  p.wgt_last() = p.wgt();
  p.E_last() = E;
  p.u_last() = p.u();
  p.r_last() = p.r();

  // Reset event variables
  p.event() = TallyEvent::KILL;
  p.event_nuclide() = NUCLIDE_NONE;
  p.event_mt() = REACTION_NONE;

  p.macro_xs().total = 0.0;
  p.macro_xs().neutron.absorption = 0.0;
  p.macro_xs().neutron.fission = 0.0;
  p.macro_xs().neutron.nu_fission = 0.0;

  // Skip void material
  if (mat_idx == -1)
    return;

  Material const& m = *materials[mat_idx];

  // Used for looping over thermal tables. Tryna use less registers
  unsigned char n_therm_tables = m.thermal_tables_.size();
  unsigned char therm_table_indx = 0;
  int next_sab_nuclide = -1;
  if (n_therm_tables) {
    next_sab_nuclide = m.thermal_tables_[therm_table_indx].index_nuclide;
  }

  unsigned i_log_union = std::log(E / energy_min_neutron) / log_spacing;

  // Add contribution from each nuclide in material
  auto const n_nuclides = m.nuclide_.size();
  for (int i = 0; i < n_nuclides; ++i) {

    auto const& i_nuclide =
      m.nuclide_[i]; // TODO test if making not a reference better
    auto* __restrict__ micro_ref {&p.neutron_xs(i_nuclide)};
    NuclideMicroXS micro;

    micro.index_sab = C_NONE;
    micro.sab_frac = 0.0;
    if (i == next_sab_nuclide) {
      const auto sab {m.thermal_tables_[therm_table_indx]};
      micro.index_sab = sab.index_table;
      micro.sab_frac = sab.fraction;
      if (E > gpu::thermal_scatt[micro.index_sab]->energy_max_)
        micro.index_sab = C_NONE;
      ++therm_table_indx;
      if (therm_table_indx < n_therm_tables) {
        next_sab_nuclide = m.thermal_tables_[therm_table_indx].index_nuclide;
      } else {
        next_sab_nuclide = -1; // done with S(a, b)
      }
    }

    NuclideMicroXS* __restrict__ use_micro = micro_ref;
    if (E != micro_ref->last_E || p.sqrtkT() != micro_ref->last_sqrtkT ||
        micro.index_sab != micro_ref->index_sab ||
        micro.sab_frac != micro_ref->sab_frac) {
      use_micro = &micro; // ummm
      auto const& nuclide = *nuclides[i_nuclide];
      micro.elastic = CACHE_INVALID;
      micro.thermal = 0.0;
      micro.thermal_elastic = 0.0;
      micro.use_ptable = false;
      micro.last_E = E;
      micro.last_sqrtkT = p.sqrtkT();

      if (nuclide.multipole_ && (E >= nuclide.multipole_->E_min_ && E <= nuclide.multipole_->E_max_)) {
        const auto& mp = *nuclide.multipole_;
        constexpr double gSQRT_PI = 1.7724538509055159927;

        double sig_s;
        double sig_a;
        double sig_f;
        
        // calculate multipole stuff...
        const double sqrtE = std::sqrt(E);
        const double invE = 1.0 / E;
        const unsigned i_window =
          std::min(static_cast<unsigned>(mp.window_info_.size() - 1),
            static_cast<unsigned>(
              (sqrtE - std::sqrt(mp.E_min_)) * mp.inv_spacing_));
        const auto& window {mp.window_info_[i_window]};

        if (p.sqrtkT() > 0.0 && window.broaden_poly) {
          // Broaden the curvefit.
          double dopp = mp.sqrt_awr_ / p.sqrtkT();
          array<double, WindowedMultipole::MAX_POLY_COEFFICIENTS>
            broadened_polynomials;

          // Broaden WMP polynomials (TODO replace with recursive version)
          const double beta = sqrtE * dopp;
          const double half_inv_dopp2 = 0.5 / (dopp * dopp);
          const double quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;
          double erf_beta;
          double exp_m_beta2;

          if (beta > 6.0) {
            // Save time, ERF(6) is 1 to machine precision.
            // beta/sqrtpi*exp(-beta**2) is also approximately 1 machine
            // epsilon.
            erf_beta = 1.;
            exp_m_beta2 = 0.;
          } else {
            erf_beta = std::erf(beta);
            exp_m_beta2 = std::exp(-beta * beta);
          }

          broadened_polynomials[0] = erf_beta / E;
          broadened_polynomials[1] = 1. / sqrtE;
          broadened_polynomials[2] =
            broadened_polynomials[0] * (half_inv_dopp2 + E) +
            exp_m_beta2 / (beta * gSQRT_PI);
          broadened_polynomials[3] =
            broadened_polynomials[1] * (E + 3.0 * half_inv_dopp2);
          const int n = mp.fit_order_ + 1;
          for (int i = 1; i < n - 3; i++) {
            double ip1_dbl = i + 1;
            broadened_polynomials[i + 3] =
              -broadened_polynomials[i - 1] * (ip1_dbl - 1.) * ip1_dbl *
                quarter_inv_dopp4 +
              broadened_polynomials[i + 1] *
                (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
          }

          for (int i_poly = 0; i_poly < mp.fit_order_ + 1; ++i_poly) {
            sig_s += mp.curvefit_(i_window, i_poly).fit_s *
                     broadened_polynomials[i_poly];
            sig_a += mp.curvefit_(i_window, i_poly).fit_a *
                     broadened_polynomials[i_poly];
            if (mp.fissionable_) {
              sig_f += mp.curvefit_(i_window, i_poly).fit_f *
                       broadened_polynomials[i_poly];
            }
          }
        } else {
          // Evaluate as if it were a polynomial
          double temp = invE;
          for (int i_poly = 0; i_poly < mp.fit_order_ + 1; ++i_poly) {
            sig_s += mp.curvefit_(i_window, i_poly).fit_s * temp;
            sig_a += mp.curvefit_(i_window, i_poly).fit_a * temp;
            if (mp.fissionable_) {
              sig_f += mp.curvefit_(i_window, i_poly).fit_f * temp;
            }
            temp *= sqrtE;
          }
        }

        // Add in pole contributions
        if (p.sqrtkT() == 0.0) {
          for (int i_pole = window.index_start; i_pole <= window.index_end;
               ++i_pole) {
            const thrust::complex<double> minus_i(0.0, -1.0);
            const thrust::complex<double> psi_chi =
              minus_i / (mp.data_[i_pole].ea - sqrtE);
            const thrust::complex<double> c_temp = psi_chi * invE;
            sig_s += (mp.data_[i_pole].rs * c_temp).real();
            sig_a += (mp.data_[i_pole].ra * c_temp).real();
            if (mp.fissionable_) {
              sig_f += (mp.data_[i_pole].rf * c_temp).real();
            }
          }
        } else {
          const double dopp = mp.sqrt_awr_ / p.sqrtkT();
          for (int i_pole = window.index_start; i_pole <= window.index_end;
               ++i_pole) {
            const thrust::complex<double> z =
              (sqrtE - mp.data_[i_pole].ea) * dopp;
            const thrust::complex<double> w_val =
              zpf8h_faddeeva(z) * dopp * invE * gSQRT_PI;
            sig_s += (mp.data_[i_pole].rs * w_val).real();
            sig_a += (mp.data_[i_pole].ra * w_val).real();
            if (mp.fissionable_) {
              sig_f += (mp.data_[i_pole].rf * w_val).real();
            }
          }
        }

        micro.total = sig_s + sig_a;
        micro.elastic = sig_s;
        micro.absorption = sig_a;
        micro.fission = sig_f;
        micro.nu_fission =
          nuclide.fissionable_
            ? micro.fission * nuclide.nu(E, EmissionMode::total)
            : 0.0;
      } else { // lookup pointwise XS

        // Find the appropriate temperature index. why would someone use
        // nearest?
        xsfloat kT = p.sqrtkT() * p.sqrtkT();
        xsfloat f;
        int i_temp;

        switch (gpu::temperature_method) {
        case TemperatureMethod::NEAREST: {
          double max_diff = INFTY;
          for (int t = 0; t < nuclide.kTs_.size(); ++t) {
            double diff = std::abs(nuclide.kTs_[t] - kT);
            if (diff < max_diff) {
              i_temp = t;
              max_diff = diff;
            }
          }
        } break;

        case TemperatureMethod::INTERPOLATION:
          // Find temperatures that bound the actual temperature
          for (i_temp = 0; i_temp < nuclide.kTs_.size() - 1; ++i_temp) {
            if (nuclide.kTs_[i_temp] <= kT && kT < nuclide.kTs_[i_temp + 1])
              break;
          }

          // Randomly sample between temperature i and i+1
          f = (kT - nuclide.kTs_[i_temp]) /
              (nuclide.kTs_[i_temp + 1] - nuclide.kTs_[i_temp]);
          if (f > prn(p.current_seed()))
            ++i_temp;
          break;
        }

        const auto& grid {nuclide.grid_[i_temp]};
        // Determine bounding indices based on which equal log-spaced
        // interval the energy is in
        int i_low = __ldg(&grid.grid_index[i_log_union]);
        int i_high = __ldg(&grid.grid_index[i_log_union + 1]) + 1;

        // Perform binary search over reduced range
        int i_grid = i_low + lower_bound_index_linear(
                               &grid.energy[i_low], &grid.energy[i_high], E);
        const auto xs_left {nuclide.xs_[i_temp][i_grid]};
        const auto xs_right {nuclide.xs_[i_temp][i_grid + 1]};
        // check for rare case where two energy points are the same
        if (grid.energy[i_grid] == grid.energy[i_grid + 1])
          ++i_grid;

        // calculate interpolation factor
        f = (E - grid.energy[i_grid]) /
            (grid.energy[i_grid + 1] - grid.energy[i_grid]);

        micro.index_temp = i_temp;
        micro.index_grid = i_grid;
        micro.interp_factor = f;

        // Calculate all microscopic cross sections
        micro.total = (1.0 - f) * xs_left.total + f * xs_right.total;
        micro.absorption =
          (1.0 - f) * xs_left.absorption + f * xs_right.absorption;

        if (nuclide.fissionable_) {
          // Calculate microscopic nuclide total cross section
          micro.fission = (1.0 - f) * xs_left.fission + f * xs_right.fission;

          // Calculate microscopic nuclide nu-fission cross section
          micro.nu_fission =
            (1.0 - f) * xs_left.nu_fission + f * xs_right.nu_fission;
        } else {
          micro.fission = 0.0;
          micro.nu_fission = 0.0;
        }

        // Calculate microscopic nuclide photon production cross section
        micro.photon_prod =
          (1.0 - f) * xs_left.photon_production + f * xs_right.photon_production;

        // Additionally calculate S(a, b) cross section data
        if (micro.index_sab >= 0) {
          int i_temp;
          xsfloat inelastic;
          xsfloat elastic;
          gpu::thermal_scatt[micro.index_sab]->calculate_xs(E, micro.last_sqrtkT,
            &i_temp, &elastic, &inelastic, p.current_seed());
          micro.thermal = micro.sab_frac * (elastic + inelastic);
          micro.thermal_elastic = micro.sab_frac * elastic;

          // calculate_elastic_xs
          if (micro.index_temp >= 0) {
            const auto& xs = nuclide.reactions_[0]->xs_[micro.index_temp].value;
            micro.elastic = (1.0 - micro.interp_factor) * xs[micro.index_grid] +
                            micro.interp_factor * xs[micro.index_grid + 1];
          }

          micro.total =
            micro.total + micro.thermal - micro.sab_frac * micro.elastic;
          micro.elastic = micro.thermal + (1.0 - micro.sab_frac) * micro.elastic;
          micro.index_temp_sab = i_temp;
        }

        // Calculate URR cross sections if needed
        if (gpu::urr_ptables_on && nuclide.urr_present_) {
          if (nuclide.urr_data_[i_temp].energy_in_bounds(E)) {
            micro.use_ptable = true;
            // TODO check storing by value
            const auto& urr = nuclide.urr_data_[i_temp];

            int i_energy = 0;
            while (E >= urr.energy_[i_energy + 1]) {
              ++i_energy;
            };

            p.stream() = STREAM_URR_PTABLE;
            double r =
              future_prn(static_cast<int64_t>(nuclide.index_), *p.current_seed());
            p.stream() = STREAM_TRACKING;

            int i_low = 0;
            while (urr.cdf_values_(i_energy, i_low) <= r) {
              ++i_low;
            };

            int i_up = 0;
            while (urr.cdf_values_(i_energy + 1, i_up) <= r) {
              ++i_up;
            };

            // Determine elastic, fission, and capture cross sections from the
            // probability table
            xsfloat elastic = 0.;
            xsfloat fission = 0.;
            xsfloat capture = 0.;
            xsfloat f;
            if (urr.interp_ == Interpolation::lin_lin) {
              // Determine the interpolation factor on the table
              f = (E - urr.energy_[i_energy]) /
                  (urr.energy_[i_energy + 1] - urr.energy_[i_energy]);

              elastic = (1. - f) * urr.xs_values_(i_energy, i_low).elastic +
                        f * urr.xs_values_(i_energy + 1, i_up).elastic;
              fission = (1. - f) * urr.xs_values_(i_energy, i_low).fission +
                        f * urr.xs_values_(i_energy + 1, i_up).fission;
              capture = (1. - f) * urr.xs_values_(i_energy, i_low).n_gamma +
                        f * urr.xs_values_(i_energy + 1, i_up).n_gamma;
            } else if (urr.interp_ == Interpolation::log_log) {
              // Determine interpolation factor on the table
              f = std::log(E / urr.energy_[i_energy]) /
                  std::log(urr.energy_[i_energy + 1] / urr.energy_[i_energy]);

              // Calculate the elastic cross section/factor
              if ((urr.xs_values_(i_energy, i_low).elastic > 0.) &&
                  (urr.xs_values_(i_energy + 1, i_up).elastic > 0.)) {
                elastic = std::exp(
                  (1. - f) * std::log(urr.xs_values_(i_energy, i_low).elastic) +
                  f * std::log(urr.xs_values_(i_energy + 1, i_up).elastic));
              } else {
                elastic = 0.;
              }

              // Calculate the fission cross section/factor
              if ((urr.xs_values_(i_energy, i_low).fission > 0.) &&
                  (urr.xs_values_(i_energy + 1, i_up).fission > 0.)) {
                fission = std::exp(
                  (1. - f) * std::log(urr.xs_values_(i_energy, i_low).fission) +
                  f * std::log(urr.xs_values_(i_energy + 1, i_up).fission));
              } else {
                fission = 0.;
              }

              // Calculate the capture cross section/factor
              if ((urr.xs_values_(i_energy, i_low).n_gamma > 0.) &&
                  (urr.xs_values_(i_energy + 1, i_up).n_gamma > 0.)) {
                capture = std::exp(
                  (1. - f) * std::log(urr.xs_values_(i_energy, i_low).n_gamma) +
                  f * std::log(urr.xs_values_(i_energy + 1, i_up).n_gamma));
              } else {
                capture = 0.;
              }
            }

            // Determine the treatment of inelastic scattering
            xsfloat inelastic = 0.;
            if (urr.inelastic_flag_ != C_NONE) {
              // get interpolation factor
              f = micro.interp_factor;

              // Determine inelastic scattering cross section
              Reaction* rx = nuclide.reactions_[nuclide.urr_inelastic_].get();
              int xs_index = micro.index_grid - rx->xs_[i_temp].threshold;
              if (xs_index >= 0) {
                inelastic = (1. - f) * rx->xs_[i_temp].value[xs_index] +
                            f * rx->xs_[i_temp].value[xs_index + 1];
              }
            }

            // Multiply by smooth cross-section if needed
            if (urr.multiply_smooth_) {
              const auto& xs = nuclide.reactions_[0]->xs_[i_temp].value;
              f = micro.interp_factor;
              micro.elastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];
              elastic *= micro.elastic;
              capture *= (micro.absorption - micro.fission);
              fission *= micro.fission;
            }

            // Check for negative values
            if (elastic < 0.) {
              elastic = 0.;
            }
            if (fission < 0.) {
              fission = 0.;
            }
            if (capture < 0.) {
              capture = 0.;
            }

            // Set elastic, absorption, fission, total, and capture x/s. Note that
            // the total x/s is calculated as a sum of partials instead of the
            // table-provided value
            micro.elastic = elastic;
            micro.absorption = capture + fission;
            micro.fission = fission;
            micro.total = elastic + inelastic + capture + fission;

            // if (simulation::need_depletion_rx) {
            //   micro.reaction[0] = capture;
            // }

            // Determine nu-fission cross-section
            if (nuclide.fissionable_) {
              micro.nu_fission =
                nuclide.nu(E, EmissionMode::total) * micro.fission;
            }
          }
        }
      }
    }

    double const& atom_density = m.atom_density_[i];
    p.macro_xs().total += atom_density * use_micro->total;
    p.macro_xs().neutron.absorption += atom_density * use_micro->absorption;
    p.macro_xs().neutron.fission += atom_density * use_micro->fission;
    p.macro_xs().neutron.nu_fission += atom_density * use_micro->nu_fission;
    if (use_micro != micro_ref) {
      *micro_ref = micro; // save stack variable back to global memory
    }
  }
}

} // namespace gpu
} // namespace openmc
