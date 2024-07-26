//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/tensor.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! UrrData contains probability tables for the unresolved resonance range.
//==============================================================================


class ContinuousURRData {
public:
  ContinuousURRData(const std::string& filename, gsl::index index);

  bool energy_in_bounds(double E) const
  {
    return energy_.front() < E && E < energy_.back();
  }

  // This takes the actual value of energy as the first argument, the temperature
  // index as the second, and the URR stream seed pointer as third. The temperature
  // is passed as an index rather than a value because the temperature grid is shared
  // across all nuclides.
  void sample(double E, int i_T, uint64_t* seed, NuclideMicroXS& xs);

private:
  vector<double> energy_; //!< incident energies
  bool has_fission_ {false};

  tensor<double, 2> alpha;
  tensor<double, 2> beta;
  tensor<double, 2> mu;
  tensor<double, 2> delta2;

  // Conditional partial values
  tensor<double, 2> nodes;
  tensor<double, 2> weights;
  tensor<double, 3> abs_values;
  tensor<double, 3> fiss_values;

  // Copy of the nuclide index for LCG stream reasons
  gsl::index index_;
};

} // namespace openmc

#endif // OPENMC_URR_H
