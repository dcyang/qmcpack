#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

add_library(minimal_hamiltonian_pool MinimalHamiltonianPool.cpp)
target_include_directories(minimal_hamiltonian_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_hamiltonian_pool PUBLIC qmcham)

set(SRC_DIR hamiltonian)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(COULOMB_SRCS test_coulomb_pbcAB.cpp test_coulomb_pbcAB_ewald.cpp test_coulomb_pbcAA.cpp
                 test_coulomb_pbcAA_ewald.cpp test_EwaldRef.cpp test_NonLocalECPotential.cpp)
set(EWALD2D_SRCS test_ewald2d.cpp test_ewald_quasi2d.cpp)
set(FORCE_SRCS test_force.cpp test_force_ewald.cpp test_stress.cpp test_spacewarp.cpp)
set(HAM_SRCS
    test_bare_kinetic.cpp
    test_density_estimator.cpp
    test_NonLocalTOperator.cpp
    test_ecp.cpp
    test_hamiltonian_pool.cpp
    test_hamiltonian_factory.cpp
    test_PairCorrEstimator.cpp
    test_SkAllEstimator.cpp
    test_QMCHamiltonian.cpp
    test_ObservableHelper.cpp
    test_Listener.cpp)

if(NOT QMC_COMPLEX)
  set(HAM_SRCS ${HAM_SRCS} test_RotatedSPOs_NLPP.cpp)
else()
  set(HAM_SRCS ${HAM_SRCS} test_SOECPotential.cpp)
endif()

set(FORCE_SRCS ${FORCE_SRCS} test_ion_derivs.cpp)

set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)

set(UTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp/pwscf.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT2} ${UTEST_DIR}/hcpBe.pwscf.h5)

set(UTEST_HDF_INPUT5 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third.h5)
maybe_symlink(${UTEST_HDF_INPUT5} ${UTEST_DIR}/C_diamond-twist-third.h5)

set(UTEST_HDF_INPUT5 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third-cart.h5)
maybe_symlink(${UTEST_HDF_INPUT5} ${UTEST_DIR}/C_diamond-twist-third-cart.h5)

foreach(fname Na2.structure.xml simple.txt)
  maybe_symlink(${CMAKE_CURRENT_SOURCE_DIR}/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(fname cn.wfnoj.xml cn.wfj.xml cn.msd-wfnoj.xml cn.msd-1det-wfnoj.xml cn.msd-1det-wfj.xml cn.msd-wfj.xml C_diamond-twist-third-cart.wfj.xml)
  maybe_symlink(${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(fname C.BFD.xml Na.BFD.xml so_ecp_test.xml C.ccECP.xml N.ccECP.xml Be.BFD.xml)
  maybe_symlink(${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(CATEGORY coulomb force ham ewald2d)
  set(UTEST_EXE test_${SRC_DIR}_${CATEGORY})
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  string(TOUPPER "${CATEGORY}_SRCS" SOURCE_FILE_VAR_NAME)
  add_executable(${UTEST_EXE} ${${SOURCE_FILE_VAR_NAME}})

  target_link_libraries(${UTEST_EXE} catch_main qmcham utilities_for_test minimal_test_pools)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} qmcwfs qmcwfs_omptarget qmcparticle qmcparticle_omptarget spline2 spline2_omptarget qmcutil platform_omptarget_LA)
  endif()

  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
endforeach()
