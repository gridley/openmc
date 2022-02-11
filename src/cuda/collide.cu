#include "openmc/cuda/atomic_agg_inc.h"
#include "openmc/cuda/collide.h"
#include "openmc/geometry.h"

namespace openmc {
namespace gpu {

__global__ void process_collision_events_device(
  EventQueueItem* __restrict__ queue,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  const unsigned p_idx = queue[tid].idx;
  Particle p(p_idx);

  p.event_collide();

  if (!p.alive()) {
    // NOTE there may be some differences with how we define CELL_BORN
    // here compared to CPU OpenMC.
    p.event_revive_from_secondary();
    if (!exhaustive_find_cell(p))
      __trap();
    if (gpu::c_micro_xs_caching)
      p.invalidate_neutron_xs();
  }

  // TODO: we should not increment this in revive_from_secondary!
  // Replace with revival from secondaries eventually
  // p.n_event()++;
  if (p.alive() && (p.material() == MATERIAL_VOID ||
                     !gpu::materials[p.material()]->fissionable_)) {
    calculate_nonfuel_xs_queue[atomicAggInc(
      &managed_calculate_nonfuel_queue_index)] = {p, p_idx};
  } else if (p.alive()) {
    calculate_fuel_xs_queue[atomicAggInc(
      &managed_calculate_fuel_queue_index)] = {p, p_idx};
  }
}

} // namespace gpu
} // namespace openmc
