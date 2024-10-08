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

add_subdirectory(CPU)

if(ENABLE_CUDA)
  add_subdirectory(CUDA)
endif()

if(ENABLE_SYCL)
  add_subdirectory(SYCL)
endif()

add_subdirectory(OMPTarget)

set(UTEST_EXE test_platforms)
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

add_executable(${UTEST_EXE} test_PlatformSelector.cpp test_AccelBLAS.cpp)
target_link_libraries(${UTEST_EXE} platform_LA platform_runtime catch_main containers)

if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} platform_omptarget_LA)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
