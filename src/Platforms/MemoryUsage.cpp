//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MemoryUsage.h"
#include <cstring>
#include <string>
#include <iomanip>
#include "Host/sysutil.h"
#include "OMPTarget/OMPallocator.hpp"
#include "MemManageAlias.hpp"

namespace qmcplusplus
{
void print_mem(const std::string& title, std::ostream& log)
{
  std::string line_separator;
  for (std::string::size_type i = 0; i < title.size() + 30; i++)
    line_separator += "=";
  log << line_separator << std::endl;
  log << "--- Memory usage report : " << title << " ---" << std::endl;
  log << line_separator << std::endl;
  log << std::right;
  log << "Available memory on node 0, free + buffers : " << std::setw(7) << (freemem() >> 20) << " MiB" << std::endl;
  log << "Memory footprint by rank 0 on node 0       : " << std::setw(7) << (memusage() >> 10) << " MiB" << std::endl;
#ifdef ENABLE_OFFLOAD
  log << "Device memory allocated via OpenMP offload : " << std::setw(7) << (getOMPdeviceMemAllocated() >> 20) << " MiB"
      << std::endl;
#endif
#ifdef ENABLE_CUDA
  log << "Device memory allocated via CUDA allocator : " << std::setw(7)
      << (compute::MemManage<PlatformKind::CUDA>::getDeviceMemAllocated() >> 20) << " MiB" << std::endl;
  log << "Free memory on the default device          : " << std::setw(7)
      << (compute::MemManage<PlatformKind::CUDA>::getDeviceFreeMem() >> 20) << " MiB" << std::endl;
#endif
#ifdef ENABLE_SYCL
  log << "Device memory allocated via SYCL allocator : " << std::setw(7)
      << (compute::MemManage<PlatformKind::SYCL>::getDeviceMemAllocated() >> 20) << " MiB" << std::endl;
  log << "Free memory on the default device          : " << std::setw(7)
      << (compute::MemManage<PlatformKind::SYCL>::getDeviceFreeMem() >> 20) << " MiB" << std::endl;
#endif
  log << line_separator << std::endl;
}

} // namespace qmcplusplus
