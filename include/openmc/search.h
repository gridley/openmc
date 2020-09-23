//! \file search.h
//! Search algorithms

#ifndef OPENMC_SEARCH_H
#define OPENMC_SEARCH_H

#include <algorithm> // for lower_bound, upper_bound
#include <ieee754.h>
#include <cmath>

namespace openmc {

// Quickly approximately calculates the natural logarithm in a portable way
inline double quicklog(double x)
{
  ieee754_double input_ieee; // Access to IEEE754 data fields of input
  ieee754_double mantissa_ieee; // Extract mantissa of the input into this variable
  input_ieee.d = x;
  mantissa_ieee.ieee.exponent = IEEE754_DOUBLE_BIAS;
  mantissa_ieee.ieee.mantissa0 = input_ieee.ieee.mantissa0;
  mantissa_ieee.ieee.mantissa1 = input_ieee.ieee.mantissa1;

  // These are a linear least-square fit of the logarithm on [1, 2], where the
  // mantissa of the double is.
  constexpr double slope = 0.682233833280656;
  constexpr double intercept = -0.637056388801094;
  constexpr double log2 = std::log(2.0);
  return slope * mantissa_ieee.d + intercept + log2 * (input_ieee.ieee.exponent - IEEE754_DOUBLE_BIAS);
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type
lower_bound_index(It first, It last, const T& value)
{
  if (*first == value) return 0;
  It index = std::lower_bound(first, last, value) - 1;
  return (index == last) ? -1 : index - first;
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type
lower_bound_index_linear(It first, It last, const T& value)
{
  if (value > *(last-1)) return last-1-first;
  It orig_last = last;
  It orig_first = first;
  while (last != first + 2) {
    unsigned interpolation = static_cast<unsigned>((value - *first) / (*(last - 1) - *first) * (last - first - 1));
    if (!interpolation) interpolation++;
    else if (interpolation == last-first) interpolation--;
    It mid = first + interpolation;
    if (*mid == value) return mid - orig_first - 1;
    else if (*mid < value) first = mid;
    else last = mid + 1;
  }
  return first - orig_first;
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type
lower_bound_index_log(It first, It last, const T& value)
{
  if (value > *(last-1)) return last-1-first;
  It orig_last = last;
  It orig_first = first;
  while (last != first + 2) {
    unsigned interpolation = static_cast<unsigned>(std::log(value/(*first)) / std::log(*(last - 1)/(*first)) * (last - first - 1));
    if (!interpolation) interpolation++;
    else if (interpolation == last-first) interpolation--;
    It mid = first + interpolation;
    if (*mid == value) return mid - orig_first - 1;
    else if (*mid < value) first = mid;
    else last = mid + 1;
  }
  return first - orig_first;
}

template<class It, class T>
typename std::iterator_traits<It>::difference_type
upper_bound_index(It first, It last, const T& value)
{
  It index = std::upper_bound(first, last, value) - 1;
  return (index == last) ? -1 : index - first;
}

} // namespace openmc

#endif // OPENMC_SEARCH_H
