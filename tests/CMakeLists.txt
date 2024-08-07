message(STATUS "Adding integration tests for QMCPACK")
message("Use --log-level=VERBOSE CMake option for details of which tests will be enabled.")

add_subdirectory(converter)
add_subdirectory(estimator)
add_subdirectory(models)
add_subdirectory(io)
add_subdirectory(performance)

if(QE_FOUND)
  add_subdirectory(pw2qmcpack)
endif()

if(RMG_FOUND)
  add_subdirectory(rmg)
  message("adding RMG tests")
endif()

message("Adding system tests for QMCPACK")

if(NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 16)
  message("TEST_MAX_PROCS was unset. Set to 16")
endif()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()

if(MIXED_PRECISION)
  set(SUCCESS_STATUS_MP FALSE)
else()
  set(SUCCESS_STATUS_MP TRUE)
endif()

if(ENABLE_OFFLOAD)
  set(SUCCESS_STATUS_OFFLOAD FALSE)
else()
  set(SUCCESS_STATUS_OFFLOAD TRUE)
endif()

if(BUILD_AFQMC)
  include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
  check_python_reqs("numpy;h5py;mpi4py;scipy;pyscf" afqmc_workflow add_tests)
  if(add_tests)
    message("Adding tests for AFQMC workflow.")
    add_subdirectory("afqmc/workflow")
  endif()
endif()

if(HAVE_PYSCF)
  if(NOT QMC_COMPLEX)
    include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
    check_python_reqs("numpy;h5py;PyscfToQmcpack;PyscfToQmcpack_Spline" pyscf_workflow add_tests)

    if(add_tests)
      message("Python dependencies met. Adding PySCF workflow tests")
      add_subdirectory("pyscf")
    endif()
  endif()
endif()

#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  add_subdirectory("heg/heg_14_gamma")
  add_subdirectory("heg/heg_54_J2rpa")
  add_subdirectory("heg/heg2d_4")
  add_subdirectory("molecules")
  add_subdirectory("solids")
endif()
