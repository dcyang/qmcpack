//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <sycl/sycl.hpp>
#include "SYCLDeviceManager.h"
#include "SYCLruntime.hpp"

namespace qmcplusplus
{
sycl::queue& getSYCLDefaultDeviceDefaultQueue() { return SYCLDeviceManager::getDefaultDeviceDefaultQueue(); }

sycl::queue createSYCLInOrderQueueOnDefaultDevice()
{
  return sycl::queue(getSYCLDefaultDeviceDefaultQueue().get_context(), getSYCLDefaultDeviceDefaultQueue().get_device(),
                     sycl::property::queue::in_order());
}

sycl::queue createSYCLQueueOnDefaultDevice()
{
  return sycl::queue(getSYCLDefaultDeviceDefaultQueue().get_context(), getSYCLDefaultDeviceDefaultQueue().get_device());
}
} // namespace qmcplusplus
