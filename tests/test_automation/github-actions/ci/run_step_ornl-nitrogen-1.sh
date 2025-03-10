#!/bin/bash

# Executes pipeline steps for the nitrogen-1 CI jobs

set -x
HOST_NAME=$(hostname -s)

case "$1" in 

  configure)

    echo "Use recent CMake v3.24.3"
    export PATH=$HOME/opt/cmake/3.24.3/bin:$PATH
    # Make current environment variables available to subsequent steps, ctest
    echo "PATH=$PATH" >> $GITHUB_ENV

    QMC_DATA_DIR=/scratch/ci/QMC_DATA_FULL

    if [ -d ${GITHUB_WORKSPACE}/../qmcpack-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../qmcpack-build
    fi

    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build"
    cd ${GITHUB_WORKSPACE}/.. && mkdir qmcpack-build && cd qmcpack-build

    # Build variants
    # Real or Complex configuration
    case "${GH_JOBNAME}" in
      *"Real"*)
        echo 'Configure for real build -DQMC_COMPLEX=0'
        IS_COMPLEX=0
      ;;
      *"Complex"*)
        echo 'Configure for complex build -DQMC_COMPLEX=1'
        IS_COMPLEX=1
      ;; 
    esac

    # Mixed or Non-Mixed (default, full) precision, used with GPU code
    case "${GH_JOBNAME}" in
      *"Mixed"*)
        echo 'Configure for mixed precision build -DQMC_MIXED_PRECISION=1'
        IS_MIXED_PRECISION=1
      ;; 
      *)
        IS_MIXED_PRECISION=0
      ;;
    esac
       
    # check the GPU architecture in use
    whoami
    groups
    /opt/rocm/llvm/bin/amdgpu-arch

    echo 'Configure for building CUDA2HIP with ROCM clang compilers'
    cmake -GNinja \
          -DCMAKE_C_COMPILER=/opt/rocm/llvm/bin/clang \
          -DCMAKE_CXX_COMPILER=/opt/rocm/llvm/bin/clang++ \
          -DQMC_MPI=0 \
          -DQMC_GPU=hip \
          -DQMC_GPU_ARCHS=gfx906 \
          -DQMC_COMPLEX=$IS_COMPLEX \
          -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DQMC_DATA=$QMC_DATA_DIR \
          ${GITHUB_WORKSPACE}
          
    ;;
  
  build)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja
    ;;
   
  test)
    echo "Running deterministic tests"
    rocm-smi --showdriverversion
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ctest --output-on-failure -L deterministic -j 32 --timeout 120 --repeat after-timeout:4
    ;;
    
esac
