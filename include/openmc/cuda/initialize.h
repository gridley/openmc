#pragma once

#include "openmc/event.h"
#include "openmc/particle.h"

#include "openmc/cuda/calculate_xs.h"

namespace openmc {
namespace gpu {

__global__ void process_initialize_events_device(unsigned queue_size,
  unsigned source_offset,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  Particle p(tid);

  if (tid < queue_size) {
    p.initialize_values();
    initialize_history(p, source_offset + tid + 1);

    if (p.alive() && (p.material() == MATERIAL_VOID ||
                       !gpu::materials[p.material()]->fissionable_)) {
      calculate_nonfuel_xs_queue[atomicAggInc(
        &managed_calculate_nonfuel_queue_index)] = {p, tid};
    } else if (p.alive()) {
      calculate_fuel_xs_queue[atomicAggInc(
        &managed_calculate_fuel_queue_index)] = {p, tid};
    }
  }
}

} // namespace gpu
} // namespace openmc
