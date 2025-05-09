#ifndef QMCPLUSPLUS_SCOPEDPROFILER_H
#define QMCPLUSPLUS_SCOPEDPROFILER_H

#include "config.h"
#if defined(ENABLE_CUDA) && !defined(QMC_CUDA2HIP)
#include "CUDA/CUDAruntime.hpp"
#include <cuda_profiler_api.h>
#endif
#ifdef USE_VTUNE_API
#include <ittnotify.h>
#endif
#ifdef USE_HPCTOOLKIT_API
#include <hpctoolkit.h>
#endif

namespace qmcplusplus
{

//non-nested
class ScopedProfiler
{
public:
  ScopedProfiler(bool active) : active_(active)
  {
    if (active_)
    {
#if defined(ENABLE_CUDA) && !defined(QMC_CUDA2HIP)
      cudaErrorCheck(cudaProfilerStart(), "cudaProfilerStart failed!");
#endif
#ifdef USE_VTUNE_API
      __itt_resume();
#endif
#ifdef USE_HPCTOOLKIT_API
      hpctoolkit_sampling_start();
#endif
    }
  }

  ~ScopedProfiler()
  {
    if (active_)
    {
#if defined(ENABLE_CUDA) && !defined(QMC_CUDA2HIP)
      cudaErrorCheck(cudaProfilerStop(), "cudaProfilerStop failed!");
#endif
#ifdef USE_VTUNE_API
      __itt_pause();
#endif
#ifdef USE_HPCTOOLKIT_API
      hpctoolkit_sampling_stop();
#endif
    }
  }

  bool isActive() const { return active_; }

private:
  const bool active_;
};
} // namespace qmcplusplus
#endif
