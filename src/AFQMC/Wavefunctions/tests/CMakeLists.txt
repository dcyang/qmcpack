#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Miguel A. Morales, LLNL
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR afqmc_wfn_factory)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/afqmc/Wavefunctions)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/choldump.h5)
SET(UTEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/wfn.dat)
SET(UTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/afqmc/C_uhf/choldump.h5)
SET(UTEST_WFN_INPUT2 ${qmcpack_SOURCE_DIR}/tests/afqmc/C_uhf/wfn_msd.dat)
SET(UTEST_HDF_INPUT3 ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/choldump_phmsd.h5)
SET(UTEST_WFN_INPUT3 ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/wfn_phmsd.dat)

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/afqmc.h5)
MAYBE_SYMLINK(${UTEST_WFN_INPUT} ${UTEST_DIR}/wfn.dat)
MAYBE_SYMLINK(${UTEST_HDF_INPUT2} ${UTEST_DIR}/afqmc_msd.h5)
MAYBE_SYMLINK(${UTEST_WFN_INPUT2} ${UTEST_DIR}/wfn_msd.dat)
MAYBE_SYMLINK(${UTEST_HDF_INPUT3} ${UTEST_DIR}/afqmc_phmsd.h5)
MAYBE_SYMLINK(${UTEST_WFN_INPUT3} ${UTEST_DIR}/wfn_phmsd.dat)

ADD_EXECUTABLE(${UTEST_EXE} test_wfn_factory.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} afqmc qmcutil
${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit;afqmc" WORKING_DIRECTORY ${UTEST_DIR})

