#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/error.h" // for writing messages
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/memory.h"
#include "openmc/nuclide.h"


#include <fmt/core.h>

#include <cmath>
#include <thrust/complex.h>
#include <algorithm> // for min

#define SQRT_PI 1.7724538509055159
#define SQRT_2 1.4142135623730951

namespace openmc {

//========================================================================
// WindowedeMultipole implementation
//========================================================================

HD inline thrust::complex<double> zpf8h_faddeeva(thrust::complex<double> z)
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


template<typename T>
HD T sgn(T val)
{
	  return (0.0 < val) - (val < 0.0);
}

// Approximates the complex E_1(x) exponential integral
HD thrust::complex<double> e1z_apprx(thrust::complex<double> z)
{
  constexpr double PI=3.141592653589793;
  constexpr double EL=0.5772156649015328;
  thrust::complex<double> ce1(0.0, 0.0);
  thrust::complex<double> cr(0.0, 0.0);
  double a0 = thrust::abs(z);
  double xt = -2.0 * std::abs(z.imag());

  // Note: not protecting z=0 case, which cannot happen in MARS
  // method (no pole has zero imaginary part)

  if (a0 <= 3.0 || z.real() < xt && a0 < 40.0) {
    // Power series
    ce1.real(1.0);
    ce1.imag(0.0);
    cr.real(1.0);
    cr.imag(0.0);

    // Note: this is truncated a bit early to get better speed.
    for (int k=1; k<50; k++) {
      cr = -cr * static_cast<double>(k) * z / std::pow(k + 1, 2);
      ce1 += cr;
      if (thrust::abs(cr) < thrust::abs(ce1) * 1e-15)
        break;
    }
    if (z.real() <= 0.0 && z.imag() == 0.0) {
      ce1 = -EL - thrust::log(-z) + z * ce1 -
            thrust::complex<double>(0.0, PI * sgn(z.imag()));
    } else {
      ce1 = -EL - thrust::log(z) + z * ce1;
    }

  } else {
    // continued fraction
    auto zd = 1.0 / z;
    auto zdc = zd;
    auto zc = zdc;
    for (int k=1; k<10; k++) {
      zd = 1.0 / (zd * static_cast<double>(k) + 1.0);
      zdc = (zd - 1.0) * zdc;
      zc += zdc;
      zd = 1.0/(zd * static_cast<double>(k) + z);
      zdc = (z * zd - 1.0) * zdc;
      zc += zdc;
      if (thrust::abs(zdc) <= thrust::abs(zc) * 1e-15 && k >= 20)
        break;
    }
    ce1 = thrust::exp(-z) * zc;
    if (z.real() <= 0.0 && z.imag() == 0.0)
      ce1 -= thrust::complex<double>(0.0, PI);
  }
  return ce1;
}


WindowedMultipole::WindowedMultipole(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  // Read scalar values.
  read_dataset(group, "spacing", inv_spacing_);
  inv_spacing_ = 1.0 / inv_spacing_;
  read_dataset(group, "sqrtAWR", sqrt_awr_);
  read_dataset(group, "E_min", E_min_);
  read_dataset(group, "E_max", E_max_);

  // Read the "data" array.  Use its shape to figure out the number of poles
  // and residue types in this data.
  xt::xtensor<std::complex<double>, 2> data_tmp;
  read_dataset(group, "data", data_tmp);
  int n_residues = data_tmp.shape()[1] - 1;

  // Read poles into *better* data format
  unsigned n_poles = data_tmp.shape()[0];
  data_.resize(n_poles);

  // Read poles to GPU-compatible struct
  for (int pole = 0; pole < data_tmp.shape()[0]; ++pole) {
    data_[pole].ea = data_tmp(pole, 0);
    data_[pole].rs = data_tmp(pole, 1);
    data_[pole].ra = data_tmp(pole, 2);
    data_[pole].rf = data_tmp(pole, 3);
  }

  // Check to see if this data includes fission residues.
  fissionable_ = (n_residues == 3);

  // Read the "windows" array and use its shape to figure out the number of
  // windows.
  xt::xtensor<int, 2> windows;
  read_dataset(group, "windows", windows);
  int n_windows = windows.shape()[0];
  windows -= 1; // Adjust to 0-based indices

  // Read the "broaden_poly" arrays.
  xt::xtensor<bool, 1> broaden_poly;
  read_dataset(group, "broaden_poly", broaden_poly);
  if (n_windows != broaden_poly.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }

  // Read the "curvefit" array.
  xt::xtensor<double, 3>
    curvefit_tmp; // Curve fit coefficients (window, poly order, reaction)
  read_dataset(group, "curvefit", curvefit_tmp);
  if (n_windows != curvefit_tmp.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }
  fit_order_ = curvefit_tmp.shape()[1] - 1;

  // Copy curvefit data into GPU-compatible memory
  std::vector<tensor<CurveFitData, 2>::size_type> cf_shape(2);
  cf_shape[0] = curvefit_tmp.shape()[0];
  cf_shape[1] = curvefit_tmp.shape()[1];
  curvefit_.resize(cf_shape);
  for (int window = 0; window < n_windows; ++window) {
    for (int poly_comp = 0; poly_comp < curvefit_tmp.shape()[1]; ++poly_comp) {
      curvefit_(window, poly_comp).fit_s = curvefit_tmp(window, poly_comp, 0);
      curvefit_(window, poly_comp).fit_a = curvefit_tmp(window, poly_comp, 1);
      curvefit_(window, poly_comp).fit_f = curvefit_tmp(window, poly_comp, 2);
    }
  }

  // Check the code is compiling to work with sufficiently high fit order
  if (fit_order_ + 1 > MAX_POLY_COEFFICIENTS) {
    fatal_error(fmt::format(
      "Need to compile with WindowedMultipole::MAX_POLY_COEFFICIENTS = {}",
      fit_order_ + 1));
  }

  // Move window information into a vector
  window_info_.resize(n_windows);
  for (int i = 0; i < n_windows; ++i) {
    window_info_[i].index_start = windows(i, 0);
    window_info_[i].index_end = windows(i, 1);
    window_info_[i].broaden_poly = broaden_poly[i];
  }
}

std::tuple<double, double, double>
WindowedMultipole::evaluate(double E, double sqrtkT) const
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(
    window_info_.size() - 1, static_cast<decltype(window_info_)::size_type>(
                               (sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Initialize the ouptut cross sections
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    array<double, MAX_POLY_COEFFICIENTS> broadened_polynomials;
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s +=
        curvefit_(i_window, i_poly).fit_s * broadened_polynomials[i_poly];
      sig_a +=
        curvefit_(i_window, i_poly).fit_a * broadened_polynomials[i_poly];
      if (fissionable_) {
        sig_f +=
          curvefit_(i_window, i_poly).fit_f * broadened_polynomials[i_poly];
      }
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly).fit_s * temp;
      sig_a += curvefit_(i_window, i_poly).fit_a * temp;
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly).fit_f * temp;
      }
      temp *= sqrtE;
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
      complx minus_i(0.0, -1.0);
      complx psi_chi = minus_i / (data_[i_pole].ea - sqrtE);
      complx c_temp = psi_chi * invE;
      sig_s += (data_[i_pole].rs * c_temp).real();
      sig_a += (data_[i_pole].ra * c_temp).real();
      if (fissionable_) {
        sig_f += (data_[i_pole].rf * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
      complx z = (sqrtE - data_[i_pole].ea) * dopp;
      complx w_val = zpf8h_faddeeva(z) * dopp * invE * SQRT_PI;
      sig_s += (data_[i_pole].rs * w_val).real();
      sig_a += (data_[i_pole].ra * w_val).real();
      if (fissionable_) {
        sig_f += (data_[i_pole].rf * w_val).real();
      }
    }
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

std::tuple<double, double, double>
WindowedMultipole::evaluate_deriv(double E, double sqrtkT) const
{
  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;
  double T = sqrtkT*sqrtkT / K_BOLTZMANN;

  if (sqrtkT == 0.0) {
    fatal_error("Windowed multipole temperature derivatives are not implemented"
      " for 0 Kelvin cross sections.");
  }

  // Locate us
  int i_window = (sqrtE - std::sqrt(E_min_)) * inv_spacing_;
  const auto& window {window_info_[i_window]};

  // Initialize the ouptut cross sections.
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // TODO Polynomials: Some of the curvefit polynomials Doppler broaden so
  // rigorously we should be computing the derivative of those.  But in
  // practice, those derivatives are only large at very low energy and they
  // have no effect on reactor calculations.

  // ==========================================================================
  // Add the contribution from the poles in this window.

  double dopp = sqrt_awr_ / sqrtkT;
  for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
    complx z = (sqrtE - data_[i_pole].ea) * dopp;
    complx w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
    sig_s += (data_[i_pole].rs * w_val).real();
    sig_a += (data_[i_pole].ra * w_val).real();
    if (fissionable_) {
      sig_f += (data_[i_pole].rf * w_val).real();
    }
  }
  double norm = -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  sig_s *= norm;
  sig_a *= norm;
  sig_f *= norm;

  return std::make_tuple(sig_s, sig_a, sig_f);
}

//========================================================================
// Non-member functions
//========================================================================

void check_wmp_version(hid_t file)
{
  if (attribute_exists(file, "version")) {
    std::array<int, 2> version;
    read_attribute(file, "version", version);
    if (version[0] != WMP_VERSION[0]) {
      fatal_error(fmt::format(
        "WMP data format uses version {}.{} whereas your installation of "
        "OpenMC expects version {}.x data.",
        version[0], version[1], WMP_VERSION[0]));
    }
  } else {
    fatal_error(fmt::format("WMP data does not indicate a version. Your "
      "installation of OpenMC expects version {}x data.", WMP_VERSION[0]));
  }
}

void read_multipole_data(int i_nuclide)
{
  // Look for WMP data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::wmp, nuc->name_});

  // If no WMP library for this nuclide, just return
  if (it == data::library_map.end()) return;

  // Check if WMP library exists
  int idx = it->second;
  std::string filename(data::libraries[idx].path_);

  // Display message
  write_message(6, "Reading {} WMP data from {}", nuc->name_, filename);

  // Open file and make sure version is sufficient
  hid_t file = file_open(filename, 'r');
  check_wmp_version(file);

  // Read nuclide data from HDF5
  hid_t group = open_group(file, nuc->name_.c_str());
  nuc->multipole_ = make_unique<WindowedMultipole>(group);
  close_group(group);
  file_close(file);
}

void broaden_wmp_polynomials(double E, double dopp, int n, double factors[])
{
  // Broadening of polynomials follows procedure outlined in C. Josey, P. Ducru,
  // B. Forget, and K. Smith, "Windowed multipole for cross section Doppler
  // broadening," J. Comput. Phys., 307, 715-727 (2016).
  // https://doi.org/10.1016/j.jcp.2015.08.013

  // Factors is already pre-allocated
  double sqrtE = std::sqrt(E);
  double beta = sqrtE * dopp;
  double half_inv_dopp2 = 0.5 / (dopp * dopp);
  double quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;

  double erf_beta;    // error function of beta
  double exp_m_beta2; // exp(-beta**2)
  if (beta > 6.0) {
    // Save time, ERF(6) is 1 to machine precision.
    // beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
    erf_beta = 1.;
    exp_m_beta2 = 0.;
  } else {
    erf_beta = std::erf(beta);
    exp_m_beta2 = std::exp(-beta * beta);
  }

  // Assume that, for sure, we'll use a second order (1/E, 1/V, const)
  // fit, and no less.
  factors[0] = erf_beta / E;
  factors[1] = 1. / sqrtE;
  factors[2] = factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 /
       (beta * SQRT_PI);
  if (n > 3) factors[3] = factors[1] * (E + 3.0 * half_inv_dopp2);

  // Perform recursive broadening of high order components (Eq. 16)
  for (int i = 1; i < n - 3; i++) {
    double ip1_dbl = i + 1;
    factors[i + 3] = -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl *
          quarter_inv_dopp4 + factors[i + 1] *
          (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
  }
}

/* Holds some values that are useful to cache throughout
 * computations. Redundantly calculating over the root
 * finding procedure would be inefficient.
 *
 * The compiler would ordinarily know to cache these
 * results if everything were confined to the same compilation
 * unit, but for sake of code organization, incomplete
 * Faddeeva calculation methods have been kept separate
 * from windowed multipole code.
 */
struct IncompleteFaddeevaCache {
  const thrust::complex<double> z;    // z
  const thrust::complex<double> wz;   // stores w(z)
  const thrust::complex<double> emz2; // stores exp(-z^2)
  const double emrz2;              // stores exp(-Re[z]^2)

  double emx2; // stores exp(-x^2)
  double erfx; // stores erf(x)
};




HD double expnexp(int n, double x)
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
HD std::pair<double, double> r_integral(double m, double z)
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
       printf("bad m\n");
    }

    double zmult = 0.5;
    const double m2inv = 1.0 / m2;

    double an = m2inv;
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
  constexpr double dawson_series[20] = {4.999995965257493680e+00, -8.333216000506398302e+01,
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
HD thrust::complex<double> jump_integral(
  IncompleteFaddeevaCache const& cache, double x)
{

  thrust::complex<double> result(0.0, 0.0);
  const double m = x - cache.z.real();

  if (std::abs(cache.z.real()) < 5.0) {
    const thrust::complex<double> beta_arg(0.0, cache.z.imag() / m);
    thrust::complex<double> b1 = -thrust::log(1.0 - beta_arg);

    // Two term recurrence for exponential term derivatives
    double a0 = 1.0;
    double a1 = 2.0 * cache.z.real();

    result += b1 * a0;
    thrust::complex<double> bpow = beta_arg;
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
    const thrust::complex<double> ii(0.0, 1.0);
    const double m2 = m*m;

    const thrust::complex<double> pp1 = 6.0 + m2*(6.0 + 3.0*m2) + zr*(m*(-3.0 - 3.0*m2) +
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

    // Note:thrust::exp(zi*(zi - 2.0*ii*zr)) = exp(-z^2) / exp(-zr^2)
    result =
      (pp2 + pp1 / (cache.emz2 / cache.emrz2 * thrust::pow(m - ii * zi, 5))) /
      (8. * std::pow(zr, 5));
  }
  return result;
}

/*
 * Documented in header file.
 */
HD thrust::complex<double> incomplete_faddeeva(
  IncompleteFaddeevaCache const& cache, double x)
{
  // avoids singularity
  if (x == cache.z.real())
    x += 1e-9;

  const double m = x - cache.z.real();
  const auto rint = r_integral(m, cache.z.real());
  const double real_part = m * rint.first + 0.5 * rint.second;

  thrust::complex<double> i_integral(
    real_part, M_PI * cache.emrz2 / cache.emx2 * (0.5 * (cache.erfx - sgn(m))));

  const auto ji = jump_integral(cache, x);

  i_integral += ji;
  i_integral *= cache.emz2 / cache.emrz2;

  thrust::complex<double> result = cache.emx2 * i_integral * thrust::complex<double>(0.0, 1.0) / M_PI;
  result += 0.5 * (cache.erfx + 1.0) * cache.wz;
  return result;
}

// Gives a good initial guess to invert the CDF of the relative speed
// distribution. See the paper for an explanation of how this works.
HD double rootfinding_bootstrap_guess(double xi, double apprx_0_cdf, double dcdx, double jump, IncompleteFaddeevaCache const& cache) {

  // TODO replace these with the Jacobian formula for the gap size which avoids 2 erf evals
  // Note that these values get moved around to simplify some logic down the line.
  double yjumplo =
    0.5 * (std::erf(cache.z.real() - 1.5 * cache.z.imag()) + 1.0);
  double yjumphi =
    0.5 * (std::erf(cache.z.real() + 1.5 * cache.z.imag()) + 1.0);

  if (xi <= apprx_0_cdf) {

    if (yjumphi > 0.5 && yjumplo < 0.5) {
      jump *= (0.5 - yjumplo) / (yjumphi - yjumplo);
      yjumphi = 0.5;
    } else if (yjumplo > 0.5) {
      yjumplo = 0.0;
      yjumphi = 0.0;
      jump = 0.0;
    }

    if (jump > apprx_0_cdf) jump = apprx_0_cdf;

    const double d = yjumphi - yjumplo;
    const double sout = (apprx_0_cdf - jump) / (0.5 - d);
    const double sinv = jump > 0.0 ? d / jump : 0.0;

    if (xi >= sout * yjumplo + jump) {
      const auto r = sout * yjumplo + jump;
      const auto a = (r - dcdx * (yjumphi - 0.5) - apprx_0_cdf) / std::pow(yjumphi - 0.5, 2);
      return 0.5 * (-dcdx + std::sqrt(std::pow(dcdx, 2) - 4.0 * (apprx_0_cdf - xi) * a)) / a + 0.5;
    } else if (xi > sout * yjumplo) {
      return (xi - sout * yjumplo) * sinv + yjumplo;
    } else {
      if (sout > 0.0)
        return xi / sout;
      else return 0.5 * yjumplo;
    }

  } else { // xi > apprx_0_cdf
    if (yjumplo < 0.5 && yjumphi > 0.5) {
      jump *= (yjumphi - 0.5) / (yjumphi - yjumplo);
      yjumplo = 0.5;
    } else if (yjumphi < 0.5) {
      yjumplo = 1.0;
      yjumphi = 1.0;
      jump = 0.0;
    }

    // Clip innapropriately large jumps
    if (jump > 1.0 - apprx_0_cdf) jump = 1.0 - apprx_0_cdf;
    const auto d = yjumphi - yjumplo;
    const auto sout = (1.0 - jump - apprx_0_cdf) / (0.5 - d);
    const auto sinv = jump > 0.0 ? d / jump : 0.0;
    const auto thresh1 = sout * (yjumplo - 0.5) + jump + apprx_0_cdf;
    const auto thresh2 = sout * (yjumplo - 0.5) + apprx_0_cdf;
    if (xi >= thresh1)
      return (xi - thresh1) / sout + yjumphi;
    else if (xi > sout * (yjumplo - 0.5) + apprx_0_cdf)
      return (xi - thresh2) * sinv + yjumplo;
    else {
      const auto a = (thresh2 - dcdx * (yjumplo - 0.5) - apprx_0_cdf)/std::pow(yjumplo - 0.5, 2);
      return 0.5 *
               (-dcdx +
                 std::sqrt(std::pow(dcdx, 2) - 4.0 * (apprx_0_cdf - xi) * a)) /
               a +
             0.5;
    }
  }

#ifdef __CUDA_ARCH__
  printf("end of bootstrap\n");
  __trap();
#endif
}


HD double WindowedMultipole::sample_target_relative_speed(
  const double& E, const double& kT, uint64_t* seed) const
{


  // Define some frequently used variables.
  const double sqrtE = std::sqrt(E);
  const double sqrtkT = std::sqrt(kT);
  const double invE = 1.0 / E;
  const double beta = sqrt_awr_ / sqrtkT; // eV^{-1/2}
  const double y = sqrtE * beta;

  // Locate window containing energy
  int i_window = min((size_t)(window_info_.size() - 1),
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Sample the effective pole
  thrust::complex<double> scat_residue(0.0, 0.0);
  thrust::complex<double> scat_pole(0.0, 0.0);
  double pole_metric = 0.0;
  for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
    const thrust::complex<double> z = data_[i_pole].ea * beta - y;
    const double this_pole_metric =
      thrust::abs(data_[i_pole].rs * zpf8h_faddeeva(z)) / (-std::log(prn(seed)));
    if (this_pole_metric > pole_metric) {
      pole_metric = this_pole_metric;
      scat_residue = -data_[i_pole].rs;
      scat_pole = data_[i_pole].ea;
    }
  }

  // Calculate location of the trough of the scattering resonance
  const double a = (thrust::complex<double>(0.0, 1.0) * (scat_residue - thrust::conj(scat_residue))).real();
  const double b = (thrust::complex<double>(0.0, 1.0) * (thrust::conj(scat_residue) * scat_pole -
                             scat_residue * thrust::conj(scat_pole)))
                     .real();
  const double c = (-(thrust::conj(scat_pole) + scat_pole)).real();
  const double d = (scat_pole * thrust::conj(scat_pole)).real();
  const double discrim = b * b - a * b * c + a * a * d;
  double s_opt; // sqrt(E) value at resonance dip
  if (discrim >= 0.0) {
    s_opt = (-b + std::sqrt(discrim)) / a;
  } else {
    s_opt = sqrtE; // ¯\_(ツ)_/¯
  }

  // Now calculate the window that the resonance dip lies in.
  // Getting the polynomial contribution within the dip is essential
  // to obtaining the correct target speed distribution.
  int i_window_pole = std::min((size_t)(window_info_.size() - 1),
    static_cast<size_t>((s_opt - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window_pole {window_info_[i_window_pole]};

  // Evaluate the polynomial part of the cross section in the resonance
  // dip at zero kelvin. While this code is somewhat repeated as in
  // ::evaluate(...), this also calculates the derivative w.r.t. sqrtE
  // as it goes, so this has therefore not been consolidated into one
  // private method.
  double polynomial_xs = 0.0;
  double polynomial_xs_slope = 0.0;
  double temp = 1.0 / (s_opt * s_opt * s_opt);
  for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
    polynomial_xs_slope +=
      curvefit_(i_window_pole, i_poly).fit_s * temp * (i_poly - 2);
    temp *= s_opt;
    polynomial_xs += curvefit_(i_window_pole, i_poly).fit_s * temp;
  }

  // Add in contribution from far away poles which negligibly influence
  // the scattering cross section:
  for (int i_pole = window_pole.index_start; i_pole <= window_pole.index_end;
       ++i_pole) {

    // Avoid duplicate poles. They are legion, and treacherous!
    if (thrust::abs(data_[i_pole].ea - scat_pole) < 1e-6)
      continue;

    thrust::complex<double> c_temp =
      -thrust::complex<double>(0.0, 1.0) / (data_[i_pole].ea - s_opt) / (s_opt * s_opt);
    polynomial_xs += (data_[i_pole].rs * c_temp).real();

    // With pole sampling, my claim is that this term shouldn't
    // actually be coming into play. Numerical evidence suggests
    // this is indeed the case.
    // polynomial_xs_slope +=
    //   (data_(i_pole, MP_RS) * c_temp *
    //     (1.0 / (data_(i_pole, MP_EA) - s_opt) - 1.0 / s_opt))
    //     .real();
  }

  // Shift and nondimensionalize the linearization.
  // The nondimensional variable "x" is centered on sqrtE.
  polynomial_xs += polynomial_xs_slope * (sqrtE - s_opt);
  polynomial_xs_slope /= beta;

  // This is the nondimensional pole passed to the incomplete Faddeeva function
  const thrust::complex<double> z = scat_pole * beta - y;

  // From here on out, we seek to solve the equation CDF(x) = xi
  const double xi = prn(seed);

  // Short circuit if resonance are far away and at sufficiently
  // high energy. Relative speed distribution reduces to a Gaussian.
  // Also, resonances don't have much of an influence if the imaginary
  // part is large.
  if (std::abs(z.real()) > 20.0 || std::abs(z.imag()) > 1.0 && y > 150.0) {
    // High energy constant cross section without resonances
    return normal_percentile(xi) / beta / SQRT_2 + sqrtE;
  } else if (thrust::abs(scat_residue) < 1e-8 || std::abs(z.real()) > 20.0 || std::abs(z.imag()) > 1.0) {
    // Low energy constant cross section without resonances
    return MARS_SAMPLE_CXS;
  }

  IncompleteFaddeevaCache cache = {
    z, zpf8h_faddeeva(z), thrust::exp(-z * z), std::exp(-z.real() * z.real()), 0.0, 0.0};

  // Normalizing factor on pole
  const double pole_term = (scat_residue * M_PI * beta * cache.wz).real();
  const double potential_term =
    0.5 * polynomial_xs * SQRT_PI * (1.0 + 2.0 * y * y) / (beta * beta);
  const double linear_term = polynomial_xs_slope * SQRT_PI * y / (beta * beta);
  const double C =
    pole_term + potential_term + linear_term; // overall normalizing constant

  const thrust::complex<double> e1 = cache.emz2 * e1z_apprx(-z * z);

  // Approximate value of w(z, 0). It has a simplified formula. The only
  // approximation comes from using a less precise complex E1(z) function;
  // the rest is exact. This is only used to kickstart the root finder, so
  // extreme precision is not required.
  const thrust::complex<double> apprx_wz0 = (0.5 * cache.wz + thrust::complex<double>(0.0, 0.5) / M_PI * e1);
  const double apprx_0_cdf = ((scat_residue * M_PI * beta * apprx_wz0).real() +
                               0.25 / std::pow(beta, 2) * polynomial_xs *
                                 (-4.0 * y + SQRT_PI * (1.0 + 2.0 * y * y)) +
                               0.5 / std::pow(beta, 2) * polynomial_xs_slope *
                                 (-(1.0 + y * y) + SQRT_PI * y)) /
                             C;

  // Approximate ∂ₓw(z, x) @ x=0
  const double dcdx =
    ((scat_residue * M_PI * beta * (thrust::complex<double>(0.0, 1.0) / M_PI / cache.z)).real() +
      std::pow(y, 2) / std::pow(beta, 2) * polynomial_xs) /
    C;

  // The approximate amount of probability gained at the resonance:
  const double jump = std::max(
    std::min((scat_residue * M_PI * beta * cache.emrz2).real() / C, 1.0), 0.0);

  // Get a good guess at the value of x. For a constant cross section
  // problem at sufficiently high energy, this is exact. It is not
  // exact for extreme low energy conditions, however, or when resonances
  // are present.
  double x = normal_percentile(rootfinding_bootstrap_guess(
               xi, apprx_0_cdf, SQRT_PI * dcdx, jump, cache)) /
             SQRT_2;
  x = std::max(-4.0, std::min(x, 4.0)); // trim to reasonable range

  // Apply up to three Newton-like corrections. The specific formula employed
  // is described in the C++ edition of the book Numerical Recipes. It
  // thresholds Halley's method when close to the root.
  for (int i = 0; i < 20; ++i) {
    cache.emx2 = std::exp(-x * x);
    cache.erfx = std::erf(x);

    // Evaluate the cumulative distribution function
    const double cdf =
      ((scat_residue * M_PI * beta * incomplete_faddeeva(cache, x)).real() +
        0.25 / (beta * beta) * polynomial_xs *
          (-2.0 * cache.emx2 * (x + 2.0 * y) +
            SQRT_PI * (1.0 + 2.0 * y * y) * (1.0 + cache.erfx)) +
        0.5 / (beta * beta) * polynomial_xs_slope *
          (-cache.emx2 * (1.0 + std::pow(x + y, 2)) +
            SQRT_PI * y * (1.0 + cache.erfx))) /
      C;
    if (std::abs(cdf - xi) < 1e-6)
      break;

    // Evaluate the probability density function (derivative for Newton)
    const double pdf =
      (cache.emx2 * (std::pow(beta, -2) * polynomial_xs * std::pow(x + y, 2) +
                      (scat_residue * beta / (cache.z - x) * thrust::complex<double>(0.0, 1.0)).real() +
                      cache.emx2 * std::pow(beta, -2) * polynomial_xs_slope *
                        std::pow(x + y, 2) * x)) /
      C;

    // Derivative of the PDF with respect to x.
    const double pdf2 =
      (-2.0 * x * pdf * C +
        cache.emx2 *
          (2.0 * std::pow(beta, -2) * polynomial_xs * (x + y) +
            (thrust::complex<double>(0.0, 1.0) * scat_residue * beta / thrust::pow(z - x, 2)).real())) /
      C;

    // Compute the Newton update
    double halley_factor =
      std::max(0.8, std::min(1.2, 1 - 0.5 * (cdf - xi) * pdf2 / pdf));
    double step = (cdf - xi) / pdf / halley_factor;
    step = std::max(-0.3, std::min(step, 0.3)); // clip step size to dx=0.3
    x -= step;
    x = std::max(-4.0, std::min(x, 4.0));
  }

  return x / beta + sqrtE;
}


} // namespace openmc
