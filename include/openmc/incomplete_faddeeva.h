#ifndef OPENMC_INCOMPLETE_FADDEEVA_H
#define OPENMC_INCOMPLETE_FADDEEVA_H

namespace openmc {

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
std::complex<double> incomplete_faddeeva(std::complex<double> z, double x);

} // namespace openmc
#endif // OPENMC_INCOMPLETE_FADDEEVA_H
