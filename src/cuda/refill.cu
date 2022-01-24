#include "openmc/cuda/refill.h"

namespace openmc {
namespace gpu {

__managed__ unsigned dead_particle_indices_indx;
__constant__ unsigned* dead_particle_indices;

__global__ void scan_for_dead_particles(unsigned n_particles)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  Particle p(tid);

  if (tid < n_particles) {
    if (!p.alive()) {
      p.event_death();
      auto loc = atomicInc(&dead_particle_indices_indx, n_particles);
      dead_particle_indices[loc] = tid;
    }
  }
}

__global__ void refill_dead_particle_slots(unsigned n_refilled,
  unsigned source_offset,
  EventQueueItem* __restrict__ calculate_nonfuel_xs_queue,
  EventQueueItem* __restrict__ calculate_fuel_xs_queue)
{
  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  unsigned p_idx = tid < n_refilled ? dead_particle_indices[tid] : 0;
  Particle p(p_idx);

  if (tid < n_refilled) {
    p.initialize_values();
    initialize_history(p, source_offset + tid + 1);

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
