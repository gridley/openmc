#include "openmc/cuda/atomic_agg_inc.h"
#include "openmc/cuda/cross_surface.h"

namespace openmc {
namespace gpu {

__global__ void process_surface_crossing_events_device(
  unsigned* __restrict__ queue, unsigned queue_size,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  const unsigned p_idx = tid < queue_size ? queue[tid] : 0;
  Particle p(p_idx);

  if (tid < queue_size) {
    p.event_cross_surface();

    // Replace with revival from secondaries eventually
    p.n_event()++;

    // These are used as booleans here, but are converted to indices shortly.
    if (p.alive() && (p.material() == MATERIAL_VOID ||
                       !gpu::materials[p.material()]->fissionable_)) {
      calculate_nonfuel_xs_queue[atomicAggInc(
        &managed_calculate_nonfuel_queue_index)] = {p, p_idx};
    } else if (p.alive()) {
      calculate_fuel_xs_queue[atomicAggInc(
        &managed_calculate_fuel_queue_index)] = {p, p_idx};
    }
  }
}

} // namespace gpu
} // namespace openmc
