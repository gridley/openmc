#pragma once

#include <cooperative_groups.h>

// from:
// https://developer.nvidia.com/blog/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/
__device__ inline unsigned atomicAggInc(unsigned* ctr)
{
  using namespace cooperative_groups;
  auto g = coalesced_threads();
  int warp_res;
  if (g.thread_rank() == 0)
    warp_res = atomicAdd(ctr, g.size());
  return g.shfl(warp_res, 0) + g.thread_rank();
}
