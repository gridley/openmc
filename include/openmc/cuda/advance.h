#pragma once

#include "openmc/cuda/atomic_agg_inc.h"
#include "openmc/cuda/calculate_xs.h" // gpu::particles
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/random_lcg.h"

namespace openmc {
namespace gpu {

extern __managed__ unsigned managed_surface_crossing_queue_index;
extern __managed__ unsigned managed_collision_queue_index;

__global__ void process_advance_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size,
  EventQueueItem* __restrict__ surface_crossing_queue,
  EventQueueItem* __restrict__ collision_queue);

} // namespace gpu
} // namespace openmc
