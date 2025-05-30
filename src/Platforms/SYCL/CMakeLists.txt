#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2022 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(SYCL_RT_SRCS SYCLDeviceManager.cpp MemManageSYCL.cpp SYCLruntime.cpp)
set(SYCL_LA_SRCS syclBLAS.cpp)

add_library(platform_sycl_runtime ${SYCL_RT_SRCS})
target_link_libraries(
  platform_sycl_runtime
  PUBLIC SYCL::host
  PRIVATE platform_host_runtime)

target_include_directories(platform_sycl_runtime PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")

add_library(platform_sycl_LA ${SYCL_LA_SRCS})
target_link_libraries(
  platform_sycl_LA
  PUBLIC platform_sycl_runtime MKL::sycl
  PRIVATE SYCL::device platform_sycl_runtime)
