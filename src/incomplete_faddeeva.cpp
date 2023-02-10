#include "openmc/incomplete_faddeeva.h"

#include <array>
#include <cmath>
#include <complex>

#include "openmc/constants.h"
#include "openmc/math_functions.h"

namespace openmc {

/* 
 *      expnexp
 *      -------
 *
 * Computes the exponential integral E_n(x) in double precision,
 * multiplied by the exponential of its argument. Treating it paired with the
 * exponential function actually saves computation.
 *
 * This is taken from the Cephes library, which doesn't have any license except
 * for a copyright. Probably just need to implement my own version of it.
 * Notably, Scipy is using this with a C adaptation, so it's probably fine to use?
 *
 * There is considerable room to improve the speed of this function. 
 * For instance, std::pow(double, double) is being called where
 * std::pow(double, int) could be called instead.
 */
double expnexp(int n, double x)
{
  double ans, r, t, yk, xk;
  double pk, pkm1, pkm2, qk, qkm1, qkm2;
  double psi, z;
  int i, k;
  constexpr double BIG = 1.44115188075855872E+17;
  constexpr double EUL = 0.57721566490153286060;
  constexpr double MACHEP = 1.111e-16;

  /* Removed Cephes error checking;
   * code calling this should never pass in a NaN, n<=0, x<=0, or n>50.
   */

  if (x > 1.0) {
    /* Continued fraction, DLMF 8.19.17 */
    k = 1;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = 1.0;
    qkm1 = x + n;
    ans = pkm1 / qkm1;

    do {
      k += 1;
      if (k & 1) {
        yk = 1.0;
        xk = n + (k - 1) / 2;
      } else {
        yk = x;
        xk = k / 2;
      }
      pk = pkm1 * yk + pkm2 * xk;
      qk = qkm1 * yk + qkm2 * xk;
      if (qk != 0) {
        r = pk / qk;
        t = fabs((ans - r) / r);
        ans = r;
      } else {
        t = 1.0;
      }
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if (fabs(pk) > BIG) {
        pkm2 /= BIG;
        pkm1 /= BIG;
        qkm2 /= BIG;
        qkm1 /= BIG;
      }
    } while (t > MACHEP);
  } else {
    /* Power series expansion, DLMF 8.19.8 */
    psi = -EUL - std::log(x);
    for (i = 1; i < n; i++) {
      psi = psi + 1.0 / i;
    }

    z = -x;
    xk = 0.0;
    yk = 1.0;
    pk = 1.0 - n;
    if (n == 1) {
      ans = 0.0;
    } else {
      ans = 1.0 / pk;
    }
    do {
      xk += 1.0;
      yk *= z / xk;
      pk += 1.0;
      if (pk != 0.0) {
        ans += yk / pk;
      }
      if (ans != 0.0)
        t = std::abs(yk / ans);
      else
        t = 1.0;
    } while (t > MACHEP);
    k = xk;

    // Note: this pow call could be an integer power instead,
    // same for the factorial here.
    ans = (std::pow(z, n - 1) * psi / std::tgamma(n)) - ans;
    ans *= std::exp(x);
  }

  return ans;
}


/*
 *      r_integral
 *      ----------
 *
 * Computes the below bivariate integral and its derivative with respect to z
 * efficiently. This is closely related to the Dawson F function, so its
 * expansion plays a pivotal role here.
 *
 *                 ∞
 *                 ⌠
 *                 ⎮    2
 *                 ⎮  -t
 *                 ⎮ e   ⋅sin(2⋅t⋅z)
 *  R(m, z)  ==    ⎮ ─────────────── dt
 *                 ⎮      2    2
 *                 ⎮     m  + t
 *                 ⌡
 *                 0
 */
std::pair<double, double> r_integral(double m, double z)
{
  double result = 0.0;
  double derivative = 0.0; // with respect to z

  const double m2 = m*m;

  /* Use asymptotic approximation for large z. The approximation
   * used for z in [-5, 5] cannot reasonably extrapolated outside.
   */
  if (std::abs(z) > 5.0) {

    // This case should never be hit in resonance upscatter calculations since only
    // the velocities over the range [-4/beta, 4/beta] are considered.
    if (std::abs(m) < 1.0) {
      std::cerr << "ERROR: Cannot call r_integral for abs(z)>5 and abs(m)<1 right now." << std::endl;
      exit(1);
    }

    double zmult = 0.5;
    const double m2inv = 1.0 / m2;

    double an = m2inv;
    double bn = 2.0;
    double cn = 2.0;
        
    // Equal to 1/z, 1/z^3, 1/z^5, ...
    double invz = 1.0 / z;

    result = zmult * an * invz;
    derivative = -result * invz;
    
    #pragma unroll
    for (int n=1; n<5; ++n) {
      an = 2.0 * n * (2.0 * n - 1.0)*an * m2inv + cn * m2inv;
      cn *= (4.0 * n + 2.0);
      zmult *= 0.25;
      invz *= 1.0/(z*z);
      const double this_term = an * zmult * invz;
      result += this_term;
      derivative -= this_term / z * (2.0 * n + 1.0);
    }
    
    return {result, derivative}; 
  }

  // Expansion of Dawson F(5 z) = c0 z + c1 z**3 + c2 z**5 + ...
  // with z in [-5, 5]. Asymptotic formula used outside.
  constexpr int n_terms = 20;
  constexpr std::array<double, n_terms> dawson_series =
     {4.999995965257493680e+00, -8.333216000506398302e+01,
      8.332309285944412522e+02, -5.948139137625151307e+03,
      3.296699363531429117e+04, -1.487301581119950279e+05,
      5.610069804904783377e+05, -1.791699162304284982e+06,
      4.855941209667988122e+06, -1.111555802822845802e+07,
      2.129946596305949986e+07, -3.378327009955811501e+07,
      4.378067747015874088e+07, -4.565152155528448522e+07,
      3.757168924686583132e+07, -2.377788318161653355e+07,
      1.113725515279613249e+07, -3.629747677761719562e+06,
      7.338263694309907733e+05, -6.922618073229778383e+04};
  const int nstar = std::min(std::max(static_cast<int>(2.72 * m2), 1), n_terms);

  // Cached value of E_nstar(m2)
  const double expn = expnexp(nstar, m2);
  double iter_expn = expn;

  const double xx = std::pow(z / 5.0, 2);

  // Backward recurse with Horner scheme
  for (int n=nstar-1; n>=1; --n) {
    iter_expn = (1.0 - n * iter_expn) / m2;
    derivative = derivative * xx + result;
    result = result * xx + iter_expn * dawson_series[n-1];
  }
  derivative *= 2.0*z*z/25.0;
  derivative += result;
  result *= z/5.0;

  // Forward recurse, iteratively build monomials
  iter_expn = expn; // reset to nstar term
  double x2 = std::pow(z / 5.0, 2*nstar-1);
  double x1 = std::pow(z / 5.0, 2*(nstar-1));
  for (int n=nstar; n<=n_terms; ++n) {
    result += dawson_series[n-1] * iter_expn * x2;
    derivative += dawson_series[n-1] * iter_expn * x1 * (2*n-1);
    x2 *= xx;
    x1 *= xx;

    // Update to next exponential integral
    iter_expn = (1.0 - m2 * iter_expn) / n;
  }

  return {result, derivative / 5.0};
}

/*
 *     jump_integral
 *     -------------
 *
 * Computes this integral efficiently. A change of variables is done such
 * that arguments to the exponential are small, and the exponential is then
 * expanded in Taylor series. This results in an infinite sum of beta
 * functions, which can be done efficiently using a recursive formula.
 *
 * The integrand becomes increasingly oscillatory as Re[z] grows, so an
 * asymptotic divergent series is employed for |Re[z]| > 5, which is at
 * most in error of about 2e-5, with rapidly diminishing error as |Re[z]|
 * grows.
 *
 *             z
 *             ⌠
 *             ⎮  ⎛ 2⎞
 *           2 ⎮  ⎝t ⎠
 *     -Re(z)  ⎮ e
 *    e        ⎮ ───── dt
 *             ⎮ x - t 
 *             ⌡
 *           Re(z)
 */
std::complex<double> jump_integral(
  IncompleteFaddeevaCache const& cache, double x)
{

  std::complex<double> result(0.0, 0.0);
  const double m = x - cache.z.real();

  if (std::abs(cache.z.real()) < 5.0) {
    const std::complex<double> beta_arg(0.0, cache.z.imag() / m);
    std::complex<double> b1 = -std::log(1.0 - beta_arg);

    // Two term recurrence for exponential term derivatives
    double a0 = 1.0;
    double a1 = 2.0 * cache.z.real();

    result += b1 * a0;
    std::complex<double> bpow = beta_arg;
    b1 -= bpow;
    result += b1 * a1 * m;
    bpow *= beta_arg;
    b1 -= 0.5 * bpow;

    // Note: this summation is highly numerically unstable
    // for |z| > 5, roughly. Hence the use of a branch to
    // handle asymptotics separately below.
    double factorial = 1.0;
    double mpow = m;
    for (double a=2.0; a<10.0; ++a) {
      double newcoe = 2.0 * cache.z.real() * a1 + 2.0 * (a - 1.0) * a0;
      a0 = a1;
      a1 = newcoe;
      factorial *= a;
      mpow *= m;
      result += b1 * newcoe * mpow / factorial;
      bpow *= beta_arg;
      b1 -= bpow / (a + 1.0);
    }
  } else { // handle asymptotic case

    // Aliases to shorten the math expression
    const double& zr = cache.z.real();
    const double& zi = cache.z.imag();
    constexpr std::complex<double> ii(0.0, 1.0);
    const double m2 = m*m;

    const std::complex<double> pp1 = 6.0 + m2*(6.0 + 3.0*m2) + zr*(m*(-3.0 - 3.0*m2) +
        zr*(m2*(2.0 + 2.0*m2) + zr*(-2.0*m2*m + 4.0*m2*m2*zr))) +
      zi*(zr*(3.0*ii + m2* (3.0*ii - 6.0*ii*m2) + zr*(m*(-4.0*ii - 4.0*ii*m2) +
              zr*(m2*(6.0*ii - 4.0*ii*m2) - 16.0*ii*m2*m*zr))) + zi*(6.0 + m2*(6.0
              - 12*m2) + zr*(m*(-3.0 - 18.0*m2) + zr*(-2 - 4.0*m2*m2 + zr*(m*(6.0
                      - 16.0*m2) - 24.0*m2*zr))) + zi*(40.0*ii*m2*m + zr*(3.0*ii
                    + m2*(18.0*ii + 4.0*ii*m2) + zr*(m*(-4.0*ii + 16.0*ii*m2) +
                      zr*(-2.0*ii + 24.0*ii*m2 + 16.0*ii*m*zr))) + zi*(3.0 + m2*(48.0 +
                      4.0*m2) + zi*(m*(-24.0*ii - 16.0*ii*m2) + zi*(-4.0 - 24.0*m2 +
                          zi*(16.0*ii*m + 4.0*zi + 4.0*ii*zr) + (-16.0*m - 4.0*zr)*zr) +
                        zr*(-24.0*ii*m2 + (-16.0*ii*m - 4.0*ii*zr)*zr)) + zr*(m*(6.0 +
                          16.0*m2) + zr*(-2.0 + 24.0*m2 + zr*(16.0*m + 4.0*zr)))))));
    const double pp2 = (-6.0 + m2*(-6.0 - 3.0*m2) + zr*(m*(3.0 + 3.0*m2) + zr*(m2*(-2.0 -
              2.0*m2) + zr*(2.0*m2*m - 4.0*m2*m2*zr))))/(m2*m2*m);

    // Note:std::exp(zi*(zi - 2.0*ii*zr)) = exp(-z^2) / exp(-zr^2)
    result =
      (pp2 + pp1 / (cache.emz2 / cache.emrz2 * std::pow(m - ii * zi, 5))) /
      (8. * std::pow(zr, 5));
  }
  return result;
}

/*
 * Documented in header file.
 */
std::complex<double> incomplete_faddeeva(
  IncompleteFaddeevaCache const& cache, double x)
{
  using namespace std::complex_literals;

  // avoids singularity
  if (x == cache.z.real())
    x += 1e-9;

  const double m = x - cache.z.real();
  const auto rint = r_integral(m, cache.z.real());
  const double real_part = m * rint.first + 0.5 * rint.second;

  std::complex<double> i_integral(
    real_part, PI * cache.emrz2 / cache.emx2 * (0.5 * (cache.erfx - sgn(m))));

  const auto ji = jump_integral(cache, x);

  i_integral += ji;
  i_integral *= cache.emz2 / cache.emrz2;

  std::complex<double> result = cache.emx2 * i_integral * 1.0i / PI;
  result += 0.5 * (cache.erfx + 1.0) * cache.wz;
  return result;
}

} // namespace openmc
