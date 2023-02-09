#ifndef OPENMC_INCOMPLETE_FADDEEVA_H
#define OPENMC_INCOMPLETE_FADDEEVA_H

#include <complex>

namespace openmc {

// Holds some values that are useful to cache throughout
// computations. Redundantly calculating over the root
// finding procedure would be inefficient.
struct IncompleteFaddeevaCache {
  const std::complex<double> z;    // z
  const std::complex<double> wz;   // stores w(z)
  const std::complex<double> emz2; // stores exp(-z^2)
  const double emrz2;              // stores exp(-Re[z]^2)

  double emx2; // stores exp(-x^2)
  double erfx; // stores erf(x)
};

/*
 * Evaluate the incomplete Faddeeva function:
 *
 *               x
 *               ⌠
 *               ⎮      2
 *           ⅈ   ⎮    -t
 *          ───  ⎮   ℯ
 *           π   ⎮  ────── dt
 *               ⎮  -t + z
 *               ⌡
 *               -∞
 */
std::complex<double> incomplete_faddeeva(
  IncompleteFaddeevaCache const& cache, double x);

} // namespace openmc
#endif // OPENMC_INCOMPLETE_FADDEEVA_H
