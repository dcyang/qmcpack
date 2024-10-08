#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2019 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//                    Bryan Clark, bclark@Princeton.edu, Princeton University
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//////////////////////////////////////////////////////////////////////////////////////

####################################
# QMC Drivers which implements various QMC algorithms
####################################

set(QMCDRIVERS
    WalkerLogBuffer.cpp
    WalkerLogCollector.cpp
    WalkerLogManager.cpp
    QMCDriverFactory.cpp
    SimpleFixedNodeBranch.cpp
    SFNBranch.cpp
    BranchIO.cpp
    MCPopulation.cpp
    QMCDriver.cpp
    QMCDriverInput.cpp
    QMCDriverNew.cpp
    WFOpt/LinearMethod.cpp
    WFOpt/ConjugateGradient.cpp
    WFOpt/Eigensolver.cpp
    WFOpt/QMCFixedSampleLinearOptimize.cpp
    WFOpt/QMCFixedSampleLinearOptimizeBatched.cpp
    WFOpt/OutputMatrix.cpp
    WFOpt/WFOptDriverInput.cpp
    WFOpt/GradientTestInput.cpp
    WFOpt/GradientTest.cpp
    Optimizers/DescentEngine.cpp
    Optimizers/HybridEngine.cpp
    WFOpt/QMCCostFunctionBase.cpp
    WFOpt/HamiltonianRef.cpp
    WFOpt/CostFunctionCrowdData.cpp
    WaveFunctionTester.cpp
    WalkerControlBase.cpp
    CloneManager.cpp
    Crowd.cpp
    QMCUpdateBase.cpp
    GreenFunctionModifiers/DriftModifierBuilder.cpp
    GreenFunctionModifiers/DriftModifierUNR.cpp
    VMC/VMCUpdatePbyP.cpp
    VMC/VMCUpdateAll.cpp
    VMC/VMCFactory.cpp
    VMC/VMCDriverInput.cpp
    VMC/SOVMCUpdatePbyP.cpp
    VMC/SOVMCUpdateAll.cpp
    DMC/DMC.cpp
    DMC/DMCUpdateAll.cpp
    DMC/DMCUpdatePbyPFast.cpp
    DMC/DMCUpdatePbyPL2.cpp
    DMC/DMCFactory.cpp
    DMC/DMCBatched.cpp
    DMC/DMCDriverInput.cpp
    DMC/DMCRefEnergy.cpp
    DMC/WalkerControlFactory.cpp
    DMC/WalkerReconfiguration.cpp
    DMC/WalkerControl.cpp
    DMC/SODMCUpdatePbyPFast.cpp
    RMC/RMC.cpp
    RMC/RMCUpdatePbyP.cpp
    RMC/RMCUpdateAll.cpp
    RMC/RMCFactory.cpp
    CorrelatedSampling/CSVMC.cpp
    CorrelatedSampling/CSVMCUpdateAll.cpp
    CorrelatedSampling/CSVMCUpdatePbyP.cpp
    CorrelatedSampling/CSUpdateBase.cpp)

set(QMCDRIVERS ${QMCDRIVERS} VMC/VMC.cpp VMC/VMCBatched.cpp WFOpt/QMCCostFunction.cpp WFOpt/QMCCostFunctionBatched.cpp)

if(HAVE_MPI)
  set(QMCDRIVERS ${QMCDRIVERS} DMC/WalkerControlMPI.cpp DMC/WalkerReconfigurationMPI.cpp)
endif(HAVE_MPI)

# if using the LMYEngine interface, make sure we add the necessary functions
if(BUILD_LMYENGINE_INTERFACE)
  set(QMCDRIVERS ${QMCDRIVERS} LMYEngineInterface/LMYE_QMCCostFunction.cpp)
  set(QMCDRIVERS ${QMCDRIVERS} LMYEngineInterface/LMYE_QMCCostFunctionBatched.cpp)
endif(BUILD_LMYENGINE_INTERFACE)

####################################
# create libqmcdriver
####################################
if(USE_OBJECT_TARGET)
  add_library(qmcdriver OBJECT ${QMCDRIVERS})
else()
  add_library(qmcdriver ${QMCDRIVERS})
endif()

target_include_directories(qmcdriver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(qmcdriver PUBLIC qmcham qmcestimators)

target_link_libraries(qmcdriver PRIVATE platform_LA Boost::boost)
if(BUILD_LMYENGINE_INTERFACE)
  target_link_libraries(qmcdriver PRIVATE formic_utils)
endif()

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
