#include "openmc/cuda/advance.h"

namespace openmc {
namespace gpu {
__managed__ unsigned managed_surface_crossing_queue_index;
__managed__ unsigned managed_collision_queue_index;

__global__ void process_advance_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size,
  EventQueueItem* __restrict__ surface_crossing_queue,
  EventQueueItem* __restrict__ collision_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  unsigned p_idx = tid < queue_size ? queue[tid].idx : 0;
  Particle p(p_idx);

  if (tid < queue_size) {
    p.event_advance();

    if (p.collision_distance() > p.boundary().distance) {
      // to surface crossing queue
      surface_crossing_queue[atomicAggInc(
        &managed_surface_crossing_queue_index)] = {p, p_idx};
    } else {
      // to collision queue
      collision_queue[atomicAggInc(&managed_collision_queue_index)] = {
        p, p_idx};
    }
  }
}

} // namespace gpu
} // namespace openmc
