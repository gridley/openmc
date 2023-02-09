#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/error.h" // for writing messages
#include "openmc/hdf5_interface.h"
#include "openmc/incomplete_faddeeva.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"
#include "openmc/random_lcg.h"

#include <fmt/core.h>

#include <algorithm> // for min
#include <cmath>

namespace openmc {

// TODO remove this
template<class UnaryOperation>
double bisection_root(double left, double right, UnaryOperation unary_op)
{

  assert(left < right);

  constexpr double TOLERANCE = 1e-6;
  constexpr double MAXITER = 100;

  // Check for differing signs
  double leftvalue = unary_op(left);
  double rightvalue = unary_op(right);

  // Sometimes it doesn't lie within -4, 4, but
  // we chop to fit in that range.
  if (leftvalue < 0.0 == rightvalue < 0.0) {
    if (rightvalue < 0.0)
      return 4.0;
    else
      return -4.0;
  }

  for (int i = 0; i < MAXITER; ++i) {

    // double middle = (right * leftvalue - left * rightvalue) / (leftvalue -
    // rightvalue);
    double middle = 0.5 * (left + right);

    double middlevalue = unary_op(middle);

    if (leftvalue < 0.0 != middlevalue < 0.0) {
      // Root in left interval
      right = middle;
      rightvalue = middlevalue;
    } else {
      // Assume to be in right interval
      left = middle;
      leftvalue = middlevalue;
    }

    if ((right - left) < TOLERANCE)
      return middle;
  }
  return (left + right) * 0.5;
}

//========================================================================
// WindowedMultipole implementation
//========================================================================

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
  read_dataset(group, "data", data_);
  int n_residues = data_.shape()[1] - 1;

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
                "array shape in WMP library for " +
                name_ + ".");
  }

  // Read the "curvefit" array.
  read_dataset(group, "curvefit", curvefit_);
  if (n_windows != curvefit_.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
                "array shape in WMP library for " +
                name_ + ".");
  }
  fit_order_ = curvefit_.shape()[1] - 1;

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

std::tuple<double, double, double> WindowedMultipole::evaluate(
  double E, double sqrtkT) const
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
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
    broaden_wmp_polynomials(
      E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s +=
        curvefit_(i_window, i_poly, FIT_S) * broadened_polynomials[i_poly];
      sig_a +=
        curvefit_(i_window, i_poly, FIT_A) * broadened_polynomials[i_poly];
      if (fissionable_) {
        sig_f +=
          curvefit_(i_window, i_poly, FIT_F) * broadened_polynomials[i_poly];
      }
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly, FIT_S) * temp;
      sig_a += curvefit_(i_window, i_poly, FIT_A) * temp;
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * temp;
      }
      temp *= sqrtE;
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = window.index_start; i_pole <= window.index_end;
         ++i_pole) {
      std::complex<double> psi_chi = -1.0i / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi * invE;
      sig_s += (data_(i_pole, MP_RS) * c_temp).real();
      sig_a += (data_(i_pole, MP_RA) * c_temp).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    for (int i_pole = window.index_start; i_pole <= window.index_end;
         ++i_pole) {
      std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
      std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
      sig_s += (data_(i_pole, MP_RS) * w_val).real();
      sig_a += (data_(i_pole, MP_RA) * w_val).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * w_val).real();
      }
    }
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

std::tuple<double, double, double> WindowedMultipole::evaluate_deriv(
  double E, double sqrtkT) const
{
  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;
  double T = sqrtkT * sqrtkT / K_BOLTZMANN;

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
    std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
    std::complex<double> w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
    sig_s += (data_(i_pole, MP_RS) * w_val).real();
    sig_a += (data_(i_pole, MP_RA) * w_val).real();
    if (fissionable_) {
      sig_f += (data_(i_pole, MP_RF) * w_val).real();
    }
  }
  double norm = -0.5 * sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  sig_s *= norm;
  sig_a *= norm;
  sig_f *= norm;

  return std::make_tuple(sig_s, sig_a, sig_f);
}

double WindowedMultipole::sample_target_relative_speed(
  const double& E, const double& kT, uint64_t* seed) const
{

  using namespace std::complex_literals;

  // Define some frequently used variables.
  const double sqrtE = std::sqrt(E);
  const double sqrtkT = std::sqrt(kT);
  const double invE = 1.0 / E;
  const double beta = sqrt_awr_ / sqrtkT; // eV^{-1/2}
  const double y = sqrtE * beta;

  // Locate window containing energy
  int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};

  // Sample the effective pole
  std::complex<double> scat_residue(0.0, 0.0);
  std::complex<double> scat_pole(0.0, 0.0);
  int selected_pole = -1;
  double pole_metric = 0.0;
  for (int i_pole = window.index_start; i_pole <= window.index_end; ++i_pole) {
    const std::complex<double> z = data_(i_pole, MP_EA) * beta - y;
    const double this_pole_metric =
      std::abs(data_(i_pole, MP_RS) / faddeeva(z));
    if (this_pole_metric > pole_metric) {
      selected_pole = i_pole;
      pole_metric = this_pole_metric;
      scat_residue = -data_(i_pole, MP_RS);
      scat_pole = data_(i_pole, MP_EA);
    }
  }

  // Calculate location of the trough of the scattering resonance
  const double a = (1.0i * (scat_residue - std::conj(scat_residue))).real();
  const double b = (1.0i * (std::conj(scat_residue) * scat_pole -
                             scat_residue * std::conj(scat_pole)))
                     .real();
  const double c = (-(std::conj(scat_pole) + scat_pole)).real();
  const double d = (scat_pole * std::conj(scat_pole)).real();
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
  int i_window_pole = std::min(window_info_.size() - 1,
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
      curvefit_(i_window_pole, i_poly, FIT_S) * temp * (i_poly - 2);
    temp *= s_opt;
    polynomial_xs += curvefit_(i_window_pole, i_poly, FIT_S) * temp;
  }

  // Add in contribution from far away poles which negligibly influence
  // the scattering cross section:
  for (int i_pole = window_pole.index_start; i_pole <= window_pole.index_end;
       ++i_pole) {

    // Avoid duplicate poles. They are legion, and treacherous!
    if (std::abs(data_(i_pole, MP_EA) - scat_pole) < 1e-6)
      continue;

    std::complex<double> c_temp =
      -1.0i / (data_(i_pole, MP_EA) - s_opt) / (s_opt * s_opt);
    polynomial_xs += (data_(i_pole, MP_RS) * c_temp).real();
    polynomial_xs_slope +=
      (data_(i_pole, MP_RS) * c_temp *
        (1.0 / (data_(i_pole, MP_EA) - s_opt) - 1.0 / s_opt))
        .real();
  }

  // Shift and nondimensionalize the linearization.
  // The nondimensional variable "x" is centered on sqrtE.
  polynomial_xs += polynomial_xs_slope * (sqrtE - s_opt);
  polynomial_xs_slope /= beta;

  // This is the nondimensional pole passed to the incomplete Faddeeva function
  const std::complex<double> z = scat_pole * beta - y;
  IncompleteFaddeevaCache cache = {
    z, faddeeva(z), std::exp(-z * z), std::exp(-z.real() * z.real()), 0.0, 0.0};

  // Normalizing factor on pole
  const double pole_term = (scat_residue * PI * beta * cache.wz).real();
  const double potential_term =
    0.5 * polynomial_xs * SQRT_PI * (1.0 + 2.0 * y * y) / (beta * beta);
  const double linear_term = polynomial_xs_slope * SQRT_PI * y / (beta * beta);
  const double C =
    pole_term + potential_term + linear_term; // overall normalizing constant

  // Print out the CDF for checking
  // for (double x=-4.0; x<4.0; x+= 0.001) {
  //   cache.emx2 = std::exp(-x * x);
  //   cache.erfx = std::erf(x);
  //   double t1 = (scat_residue * PI * beta * incomplete_faddeeva(cache,
  //   x)).real(); double t2 = 0.25 /(beta*beta) * polynomial_xs * (-2.0 *
  //   std::exp(-x*x)*(x+2.0*y)+SQRT_PI*(1.0+2.0*y*y)*(1.0+std::erf(x))); double
  //   t3 = 0.5 /(beta*beta) * polynomial_xs_slope *
  //   (-std::exp(-x*x)*(1.0+std::pow(x+y, 2))+SQRT_PI*y*(1.0+std::erf(x)));
  //   std::cout << x << "    " << (t1+t2+t3)/C << std::endl;
  // }
  // exit(0);

  // Find the nondimensional velocity with inverse CDF sampling. TODO
  // experiment with a more efficient root finding method.
  const double xi = prn(seed);
  auto cdf = [=](double x) mutable {
    cache.emx2 = std::exp(-x * x);
    cache.erfx = std::erf(x);
    return ((scat_residue * PI * beta * incomplete_faddeeva(cache, x)).real() +
             0.25 / (beta * beta) * polynomial_xs *
               (-2.0 * cache.emx2 * (x + 2.0 * y) +
                 SQRT_PI * (1.0 + 2.0 * y * y) * (1.0 + cache.erfx)) +
             0.5 / (beta * beta) * polynomial_xs_slope *
               (-cache.emx2 * (1.0 + std::pow(x + y, 2)) +
                 SQRT_PI * y * (1.0 + cache.erfx))) /
           C;
  };
  auto shifted_cdf = [xi, cdf](double x) mutable { return cdf(x) - xi; };
  const double x_sample = bisection_root(-4.0, 4.0, shifted_cdf);

  return x_sample / beta + sqrtE;
}

//========================================================================
// Non-member functions
//========================================================================

void check_wmp_version(hid_t file)
{
  if (attribute_exists(file, "version")) {
    array<int, 2> version;
    read_attribute(file, "version", version);
    if (version[0] != WMP_VERSION[0]) {
      fatal_error(fmt::format(
        "WMP data format uses version {}.{} whereas your installation of "
        "OpenMC expects version {}.x data.",
        version[0], version[1], WMP_VERSION[0]));
    }
  } else {
    fatal_error(fmt::format("WMP data does not indicate a version. Your "
                            "installation of OpenMC expects version {}x data.",
      WMP_VERSION[0]));
  }
}

void read_multipole_data(int i_nuclide)
{
  // Look for WMP data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::wmp, nuc->name_});

  // If no WMP library for this nuclide, just return
  if (it == data::library_map.end())
    return;

  // Check if WMP library exists
  int idx = it->second;
  std::string& filename = data::libraries[idx].path_;

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
  factors[2] =
    factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 / (beta * SQRT_PI);
  if (n > 3)
    factors[3] = factors[1] * (E + 3.0 * half_inv_dopp2);

  // Perform recursive broadening of high order components (Eq. 16)
  for (int i = 1; i < n - 3; i++) {
    double ip1_dbl = i + 1;
    factors[i + 3] =
      -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl * quarter_inv_dopp4 +
      factors[i + 1] * (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
  }
}

} // namespace openmc
