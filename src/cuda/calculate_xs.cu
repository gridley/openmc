#include "openmc/cuda/calculate_xs.h"

// TODO potentially load p.sqrtkT() to a register

namespace openmc {
namespace gpu {

__constant__ unique_ptr<Material>* materials;
__constant__ unique_ptr<Nuclide>* nuclides;
__constant__ unique_ptr<ThermalScattering>* thermal_scatt;
__constant__ Particle* particles;
__constant__ xsfloat energy_min_neutron;
__constant__ xsfloat energy_max_neutron;
__constant__ xsfloat log_spacing;
__constant__ unsigned number_nuclides;
__constant__ bool need_depletion_rx;

__managed__ unsigned managed_calculate_fuel_queue_index;
__managed__ unsigned managed_calculate_nonfuel_queue_index;

} // namespace gpu
} // namespace openmc
