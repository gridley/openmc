#pragma once

#include <thrust/complex.h>
#include "openmc/event.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/thermal.h" // ThermalScattering
#include "openmc/geometry.h"
#include "openmc/reaction_product.h" // EmissionMode
#include "openmc/search.h"
#include "openmc/settings.h" // BLOCKSIZE
#include "openmc/wmp.h" // TODO remove?
#include "openmc/cuda/calculate_xs.h"

namespace openmc {
namespace gpu {

__device__ double sample_nig(double alpha, double beta, double mu, double delta2, uint64_t* seed) {

  double u1 = prn(seed);
  double u2 = prn(seed);
  double R = std::sqrt(-2.0 * std::log(u1));
  double phi = 2.0 * M_PI * u2;

  // Two independent normal variates are obtained:
  double nv1 = R * std::cos(phi);
  double nv2 = R * std::sin(phi);

  // Sample the inverse Gaussian distribution, z
  double mu_ig = std::sqrt(delta2 / (alpha*alpha-beta*beta));
  double w = mu_ig * nv1 * nv1;
  double c = 0.5 * mu_ig / delta2;
  double z = mu_ig + c * (w - std::sqrt(w*(4*delta2+w)));
  if (prn(seed) >= mu_ig / (mu_ig + z)) {
    z = mu_ig * mu_ig / z;
  }

  // Sample the NIG distribution, which is a mixture over IGs
  return std::sqrt(z) * nv2 + beta * z + mu;
}

__device__ inline thrust::complex<double> zpf8h_faddeeva(thrust::complex<double> z)
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

// This class is a NuclideMicroXS (basically) when we want to use the
// UseMicroCache=true option below, and otherwise contains nothing
template<bool UseMicroCache>
struct NuclideMicroXSDummy {};

// Full micro cache on the stack
template<>
struct NuclideMicroXSDummy<true> : NuclideMicroXS {};

// Reduced data micro data on the stack (not written
// back to global memory ever)
template<>
struct NuclideMicroXSDummy<false> {
  int index_sab {-1};
  int index_temp {-1};
  int index_temp_sab {-1};
  int index_grid {-1};
  xsfloat sab_frac {0.0};
  xsfloat total {0.0};
  xsfloat elastic {0.0};
  xsfloat absorption {0.0};
  xsfloat fission {0.0};
  xsfloat nu_fission {0.0};
  xsfloat interp_factor {0.0};
  xsfloat thermal {0.0};
  xsfloat thermal_elastic {0.0};
};

// It makes the pointwise code _slightly_ faster to known whether WMP
// will not be used in advance. Probably some register optimizations being
// done by the CUDA compiler.
//
// ForCollision -- runs right before a collision and caches all micro XS,
// in cache-free mode.
template <bool UseWMP, bool UseMicroCache, bool ForCollision = false>
__global__ void  process_calculate_xs_events_device_wmp(
  EventQueueItem* __restrict__ queue)
{
  using EmissionMode = ReactionProduct::EmissionMode;
  const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  const unsigned idx = queue[tid].idx;
  const xsfloat E = __ldg(&queue[tid].E); // is ldg actually doing much for us here?
  const int mat_idx = __ldg(&queue[tid].material);
  double cutoff; // used only for pre-collision. Compiler will eliminate otherwise
  double prob; // ^^^^^
  Particle p(idx);

  static_assert(UseMicroCache == false);

  // Store pre-collision particle properties
  // TODO potentially remove this stuff???
  // In fact, do we even need this at all for the purposes of what I hope to achieve?
  if constexpr (!ForCollision) {
    p.wgt_last() = p.wgt();
    p.E_last() = E;
    p.u_last() = p.u();
    p.r_last() = p.r();

    // Reset event variables
    p.event() = TallyEvent::KILL;
    p.event_nuclide() = NUCLIDE_NONE;
    p.event_mt() = REACTION_NONE;

    // TODO potentially put this on teh stack (also advantageous
    // for implementation of cache-free collisions)
    p.macro_xs().total = 0.0;
    p.macro_xs().neutron.absorption = 0.0;
    p.macro_xs().neutron.fission = 0.0;
    p.macro_xs().neutron.nu_fission = 0.0;
  } else {
    cutoff = prn(p.current_seed()) * p.macro_xs().total;
    prob = 0.0;
    p.event_nuclide() = 0; // floating point issue protection (maybe remove?)
  }


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
  for (int i = 0; i < number_nuclides; ++i) { // loop over global nuclide array index

    // Continue through material until we reach the
    int i_nuclide = m.mat_nuclide_index_[i]; // material's nuclide index
    if constexpr (ForCollision) {
      // Check if we are done finding the collision nuclide, but still
      // need to hit the syncthreads below to avoid locking.
      if (cutoff == 1e9) i_nuclide = -1;
    }

    if (i_nuclide != -1) { // if global nuclide index present in material, enter this block

      NuclideMicroXS* __restrict__ micro_ref {nullptr};
      NuclideMicroXSDummy<UseMicroCache> micro;

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

      auto const& nuclide = *nuclides[i];

      micro.thermal = 0.0;
      micro.thermal_elastic = 0.0;
      micro.elastic = CACHE_INVALID;


      // Find the appropriate temperature index. why would someone use
      // nearest?
      xsfloat kT = p.sqrtkT() * p.sqrtkT();

      if (gpu::urr_ptables_on && nuclide.has_urr_ && nuclide.urr.energy_in_bounds(E)) {
        double T = kT / K_BOLTZMANN;
        constexpr double urr_temperatures[6] = {250.0, 294.0, 600.0, 900.0, 1200.0, 2500.0};
        int i_T;
        for (i_T=0; i_T<6; ++i_T) {
          if (urr_temperatures[i_T] > T) {
            break;
          }
        }

        int energy_index = lower_bound_index(nuclide.urr.energy_.begin(), nuclide.urr.energy_.end(), E);

	// Form a unique RNG seed that corresponds perfectly to the particle's current energy and nuclide
	double baseval = static_cast<double>(1e12 * i_nuclide) + E * 1e8;
	uint64_t * fseed = reinterpret_cast<uint64_t*>(&baseval);
	// Advance the prn once to remove any possible correlations
	prn(fseed);

	// get the parameters for the distribution
	double f = (E - nuclide.urr.energy_[energy_index]) /
	        (nuclide.urr.energy_[energy_index + 1] - nuclide.urr.energy_[energy_index]);
 
        double a = (1.0-f)*  nuclide.urr.alpha(energy_index, i_T) + f*nuclide.urr.alpha(energy_index+1, i_T);
        double b = (1.0-f)*  nuclide.urr.beta(energy_index, i_T) +  f*nuclide.urr.beta(energy_index+1, i_T);
        double m = (1.0-f)*  nuclide.urr.mu(energy_index, i_T)   +  f*nuclide.urr.mu(energy_index+1, i_T);
        double d2 = (1.0-f)* nuclide.urr.delta2(energy_index, i_T)+ f*nuclide.urr.delta2(energy_index+1, i_T);
        micro.total = sample_nig(0.5 * (a + b), -0.5 * (b - a), m, d2, fseed);

	// Compute the conditional expectations of the partials
        constexpr int bary_order = 5;
        double denom = 0.0;
        double abs = 0.0; // note "abs" is actually referring to capture here
        double fiss = 0.0;
      
        // Evaluated at next higher grid
        double abs2 = 0.0;
        double fiss2 = 0.0;
      
        for (int j=0; j<bary_order; ++j) {
          double term = nuclide.urr.weights(energy_index, j) / (micro.total - nuclide.urr.nodes(energy_index, j));
          abs += nuclide.urr.abs_values(energy_index, j, i_T) * term;
          if (nuclide.urr.has_fission_)
            fiss += nuclide.urr.fiss_values(energy_index, j, i_T) * term;
          denom += term;
        }
        abs /= denom;
        fiss /= denom;
      
        denom = 0.0;
        for (int j=0; j<bary_order; ++j) {
          double term = nuclide.urr.weights(energy_index+1, j) / (micro.total - nuclide.urr.nodes(energy_index+1, j));
          abs2 += nuclide.urr.abs_values(energy_index+1, j, i_T) * term;
          if (nuclide.urr.has_fission_)
            fiss2 += nuclide.urr.fiss_values(energy_index+1, j, i_T) * term;
          denom += term;
        }
        abs2 /= denom;
        fiss2 /= denom;
      
        micro.fission = fiss * (1.0 - f) + fiss2 * f;

        micro.absorption = abs * (1.0 - f) + abs2 * f + micro.fission;
        micro.elastic = micro.total - micro.absorption;

        // Determine nu-fission cross-section
        if (nuclide.fissionable_) {
          micro.nu_fission =
            nuclide.nu(E, EmissionMode::total) * micro.fission;
        }

        if (micro.total < 0.0 || isnan(micro.total) || isnan(abs) || isnan(fiss)) {
          micro.total = 1e-6;
          micro.absorption = 1e-6;
          micro.fission = 0.0;
          micro.nu_fission = 0.0;
          micro.elastic = 0.0;
          return;
        }

	// Just gives some reasonable values here
	//
	micro.index_temp = 0;
	micro.index_grid = 0;
      } else {

        switch (gpu::temperature_method) {
        case TemperatureMethod::NEAREST: {
          xsfloat max_diff = INFTY;
          for (int t = 0; t < nuclide.kTs_.size(); ++t) {
            xsfloat diff = std::abs(nuclide.kTs_[t] - kT);
            if (diff < max_diff) {
              micro.index_temp = t;
              max_diff = diff;
            }
          }
        } break;

        case TemperatureMethod::INTERPOLATION:
          // Find temperatures that bound the actual temperature
          for (micro.index_temp = 0; micro.index_temp < nuclide.kTs_.size() - 1; ++micro.index_temp) {
            if (nuclide.kTs_[micro.index_temp] <= kT && kT < nuclide.kTs_[micro.index_temp + 1])
              break;
          }

          // Randomly sample between temperature i and i+1
          micro.interp_factor = (kT - nuclide.kTs_[micro.index_temp]) /
              (nuclide.kTs_[micro.index_temp + 1] - nuclide.kTs_[micro.index_temp]);
          if (micro.interp_factor > prn(p.current_seed()))
            ++micro.index_temp;
          break;
        }

        const auto& grid {nuclide.grid_[micro.index_temp]};
        // Determine bounding indices based on which equal log-spaced
        // interval the energy is in
        int i_low = __ldg(&grid.grid_index[i_log_union]);
        int i_high = __ldg(&grid.grid_index[i_log_union + 1]) + 1;

        // Perform binary search over reduced range
        micro.index_grid = i_low + lower_bound_index_linear(
                               &grid.energy[i_low], &grid.energy[i_high], E);
        const auto xs_left {nuclide.xs_[micro.index_temp][micro.index_grid]};
        const auto xs_right {nuclide.xs_[micro.index_temp][micro.index_grid + 1]};
        // check for rare case where two energy points are the same
        if (grid.energy[micro.index_grid] == grid.energy[micro.index_grid + 1])
          ++micro.index_grid;

        // calculate interpolation factor
        micro.interp_factor = (E - grid.energy[micro.index_grid]) /
            (grid.energy[micro.index_grid + 1] - grid.energy[micro.index_grid]);

        // Calculate all microscopic cross sections
        micro.total = (1.0 - micro.interp_factor) * xs_left.total + micro.interp_factor * xs_right.total;
        micro.absorption =
          (1.0 - micro.interp_factor) * xs_left.absorption + micro.interp_factor * xs_right.absorption;

        if (nuclide.fissionable_) {
          // Calculate microscopic nuclide total cross section
          micro.fission = (1.0 - micro.interp_factor) * xs_left.fission + micro.interp_factor * xs_right.fission;

          // Calculate microscopic nuclide nu-fission cross section
          micro.nu_fission =
            (1.0 - micro.interp_factor) * xs_left.nu_fission + micro.interp_factor * xs_right.nu_fission;
        } else {
          micro.fission = 0.0;
          micro.nu_fission = 0.0;
        }

        // Additionally calculate S(a, b) cross section data
        xsfloat thermal = 0.0;
        if (micro.index_sab >= 0) {
          int i_temp;
          xsfloat inelastic;
          xsfloat elastic;
          // TODO cache sqrtkT earlier on the stack? Used to use micro.last_sqrtkT here.
          gpu::thermal_scatt[micro.index_sab]->calculate_xs(E, p.sqrtkT(),
            &i_temp, &elastic, &inelastic, p.current_seed());
          thermal = micro.sab_frac * (elastic + inelastic);
          micro.thermal = thermal;
          micro.thermal_elastic = micro.sab_frac * elastic;

          // calculate_elastic_xs
          if (micro.index_temp >= 0) {
            const auto& xs = nuclide.reactions_[0]->xs_[micro.index_temp].value;
            micro.elastic = (1.0 - micro.interp_factor) * xs[micro.index_grid] +
                            micro.interp_factor * xs[micro.index_grid + 1];
          }

          micro.total =
            micro.total + thermal - micro.sab_frac * micro.elastic;
          micro.elastic = thermal + (1.0 - micro.sab_frac) * micro.elastic;
          micro.index_temp_sab = i_temp;
        }
      }


        // TODO remove reference here probably
        double const& atom_density = m.atom_density_[i_nuclide];
        if constexpr (!ForCollision) {
          p.macro_xs().total += atom_density * micro.total;
          p.macro_xs().neutron.absorption += atom_density * micro.absorption;
          p.macro_xs().neutron.fission += atom_density * micro.fission;
          p.macro_xs().neutron.nu_fission += atom_density * micro.nu_fission;
        } else { // ForCollision
          prob += atom_density * micro.total;
          if (prob >= cutoff) {
            p.event_nuclide() = i; // TODO can remove this..
            // TODO make this not suck!
            NuclideMicroXS onstack;
            onstack.index_sab = micro.index_sab;
            onstack.index_temp = micro.index_temp;
            onstack.index_temp_sab = micro.index_temp_sab;
            onstack.index_grid = micro.index_grid;
            onstack.sab_frac = micro.sab_frac;
            onstack.total = micro.total;
            onstack.elastic = micro.elastic;
            onstack.absorption = micro.absorption;
            onstack.fission = micro.fission;
            onstack.nu_fission = micro.nu_fission;
            onstack.interp_factor = micro.interp_factor;
            onstack.thermal = micro.thermal;
            onstack.thermal_elastic = micro.thermal_elastic;
            p.neutron_xs(0) = onstack;
            cutoff = 1e9; // prevent any more updates
            queue[tid].material = i; // for sorting collision nuclide
          }
        }
      }
      __syncwarp();
    }

  if constexpr (ForCollision) {
    if (cutoff != 1e9) {
	    p.event_nuclide() = 0;
          NuclideMicroXS onstack;
          onstack.index_sab = 0;
          onstack.index_temp = 0;
          onstack.index_temp_sab = 0;
          onstack.index_grid = 0;
          onstack.sab_frac = 0.0;
          onstack.total = 1e-6;
          onstack.elastic = 1e-6;
          onstack.absorption = 0.0;
          onstack.fission = 0.0;
          onstack.nu_fission = 0.0;
          onstack.interp_factor = 0.0;
          onstack.thermal = 0.0;
          onstack.thermal_elastic = 0.0;
          p.neutron_xs(0) = onstack;
    }
  }
}

}
}
