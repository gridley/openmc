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
  ContinuousURRData(const std::string& filename, int index);
  ContinuousURRData() = default;
  ContinuousURRData(ContinuousURRData&& other) noexcept = default;
  ContinuousURRData& operator=(ContinuousURRData&& other) noexcept = default;

  HD bool energy_in_bounds(double E) const
  {
    return Elo_ < E && E < Ehi_;
  }

  // This is manually inlined in the XS lookup kernel
  // HD void sample(double E, int i_T, uint64_t* seed, NuclideMicroXS& xs);

  // These are cached to avoid unnecessary global memory accesses
  double Elo_ {0.0};
  double Ehi_ {0.0};

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
  int index_;
};

} // namespace openmc

#endif // OPENMC_URR_H
