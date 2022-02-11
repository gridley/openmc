#include "openmc/cuda/advance.h"
#include "openmc/geometry.h"

namespace openmc {
namespace gpu {
__managed__ unsigned managed_surface_crossing_queue_index;
__managed__ unsigned managed_collision_queue_index;

__global__ void process_advance_events_device(
  unsigned* __restrict__ queue, unsigned* __restrict__ surface_crossing_queue,
  EventQueueItem* __restrict__ collision_queue)
{
  const unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  const auto p_idx = queue[tid];
  Particle p(p_idx);

  // Find the distance to the nearest boundary
  const BoundaryInfo bdry = distance_to_boundary(p);

  const double collision_distance =
    p.macro_xs().total == 0.0 ? INFINITY : -std::log(prn(p.current_seed())) / p.macro_xs().total;
  const double distance = std::min(bdry.distance, collision_distance);

  // Advance particle
  const auto nc = p.n_coord();
  for (int j = 0; j < nc; ++j) {
    p.coord(j).r += distance * p.coord(j).u;
  }

  p.boundary() = bdry;

  // Score track-length estimate of k-eff
  // if (type() == ParticleType::neutron) {
  p.keff_tally_tracklength() += p.wgt() * distance * p.macro_xs().neutron.nu_fission;
  // }

  if (collision_distance > bdry.distance) {
    // to surface crossing queue
    surface_crossing_queue[atomicAggInc(
      &managed_surface_crossing_queue_index)] = p_idx;
  } else {
    // to collision queue
    collision_queue[atomicAggInc(&managed_collision_queue_index)] = {p, p_idx};
  }

}

} // namespace gpu
} // namespace openmc
