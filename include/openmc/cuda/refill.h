#pragma once

#include "openmc/cuda/atomic_agg_inc.h"
#include "openmc/cuda/calculate_xs.h"
#include "openmc/particle.h"
#include "openmc/simulation.h" // initialize_history
#include "openmc/event.h" // EventQueueItem

namespace openmc {
namespace gpu {

extern __managed__ unsigned dead_particle_indices_indx;
extern __constant__ unsigned* dead_particle_indices;

__global__ void scan_for_dead_particles(unsigned n_particles);

__global__ void refill_dead_particle_slots(unsigned n_refilled,
  unsigned source_offset,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue);

} // namespace gpu
} // namespace openmc
