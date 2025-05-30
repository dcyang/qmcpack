//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <memory>
#include <iostream>
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/MemManageCUDA.hpp"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{
TEST_CASE("CUDA_allocators", "[CUDA]")
{
  { // CUDAManagedAllocator
    Vector<double, CUDAManagedAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()), "cudaPointerGetAttributes failed!");
#if (CUDART_VERSION >= 10000 || HIP_VERSION_MAJOR >= 6)
    REQUIRE(attr.type == cudaMemoryTypeManaged);
#endif
  }
  { // CUDAAllocator
    Vector<double, CUDAAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()), "cudaPointerGetAttributes failed!");
#if (CUDART_VERSION >= 10000 || HIP_VERSION_MAJOR >= 6)
    REQUIRE(attr.type == cudaMemoryTypeDevice);
#else
    REQUIRE(attr.memoryType == cudaMemoryTypeDevice);
#endif
  }
  { // CUDAHostAllocator
    Vector<double, CUDAHostAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()), "cudaPointerGetAttributes failed!");
#if (CUDART_VERSION >= 10000 || HIP_VERSION_MAJOR >= 6)
    REQUIRE(attr.type == cudaMemoryTypeHost);
#else
    REQUIRE(attr.memoryType == cudaMemoryTypeHost);
#endif
  }
#if !defined(QMC_DISABLE_HIP_HOST_REGISTER)
  { // CUDALockedPageAllocator
    Vector<double, CUDALockedPageAllocator<double>> vec(1024);
    cudaPointerAttributes attr;
    cudaErrorCheck(cudaPointerGetAttributes(&attr, vec.data()), "cudaPointerGetAttributes failed!");
#if (CUDART_VERSION >= 10000 || HIP_VERSION_MAJOR >= 6)
    REQUIRE(attr.type == cudaMemoryTypeHost);
#else
    REQUIRE(attr.memoryType == cudaMemoryTypeHost);
#endif
    Vector<double, CUDALockedPageAllocator<double>> vecb(vec);
  }
#endif
  { // CUDALockedPageAllocator zero size and copy constructor
    Vector<double, CUDALockedPageAllocator<double>> vec;
    Vector<double, CUDALockedPageAllocator<double>> vecb(vec);
  }
}

} // namespace qmcplusplus
