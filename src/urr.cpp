#include "openmc/urr.h"

#include <algorithm> // any_of
#include <array>
#include <iostream>

namespace openmc {

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

// Samples an inverse Gaussian random variable
HD double sample_ig(double mu, double lam, uint64_t* seed) {

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
HD double sample_nig(double alpha, double beta, double mu, double delta2, uint64_t* seed) {
  // assert(std::abs(alpha) > std::abs(beta));
  double z = sample_ig(std::sqrt(delta2 / (alpha*alpha-beta*beta)), delta2, seed);
  return std::sqrt(z) * normal_variate(0.0, 1.0, seed) + beta * z + mu;
}

void ContinuousURRData::sample(double E, int i_T, uint64_t* seed, NuclideMicroXS& xs) {

  // TODO use linear index
  int energy_index = lower_bound_index(energy_.begin(), energy_.end(), E);

  // We use more random numbers than in the single table case. It's
  // an undefined number but quite likely less than 400.
  // TODO use a better method for sampling a normal RV
  uint64_t fseed = future_seed(static_cast<uint64_t>(400 * index_), *seed);

  // energy interpolation factor. Using stochastic interpolation.
  double f = (E - energy_[energy_index]) /
        (energy_[energy_index + 1] - energy_[energy_index]);
  if (prn(&fseed) < f) energy_index++;


  double a = alpha(energy_index, i_T);
  double b = beta(energy_index, i_T);
  double m = mu(energy_index, i_T);
  double d2 =delta2(energy_index, i_T);

  double sigt = sample_nig(0.5 * (a + b), -0.5 * (b - a), m, d2, &fseed);

  if (sigt <= 0.0 || isnan(sigt) ) {
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
  if (abs < 0.0) abs = 0.0;
    // printf("negative absorption!\n");
  if (fiss < 0.0) fiss = 0.0;
    // printf("negative fiss!\n");

  if (isnan(abs) || isnan(fiss)) {
    sigt = 1e-4;
    xs.total = 1e-4;
    xs.absorption = 1e-4;
    xs.fission = 0.0;
    xs.nu_fission = 0.0;
    xs.elastic = 0.0;
    return;
  }

  // #pragma omp critical
  //   {
  //     double fval = fiss == 0.0? 1.0 : fiss / averages(energy_index, i_T, 2);
  //     printf("%f %f %f\n", sigt / averages(energy_index, i_T, 0),
  //         abs / averages(energy_index, i_T, 1),
  //         fval);
  //   }

  // Get the inelastic, whatever other reaction contributions
  double non_abs_non_el = xs.total - xs.absorption - xs.elastic;


  // This is a direct lookup from the table. The values
  // do not necessarily have averages equal to the pointwise
  // grid given in file 3.
  xs.absorption = abs + fiss;
  xs.fission = fiss;
  xs.elastic = sigt - abs - fiss;

  xs.total = xs.absorption + xs.elastic + non_abs_non_el;
  if (isnan(xs.total)) {
    sigt = 1e-4;
    xs.total = 1e-4;
    xs.absorption = 1e-4;
    xs.fission = 0.0;
    xs.nu_fission = 0.0;
    xs.elastic = 0.0;
    return;
  }

}

} // namespace openmc
