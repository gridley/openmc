#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/calculate_xs.h"

namespace openmc {
namespace gpu {

__global__ void process_collision_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue);

} // namespace gpu
} // namespace openmc
