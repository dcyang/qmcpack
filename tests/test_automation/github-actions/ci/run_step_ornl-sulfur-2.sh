#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in

  # Configure qmcpack using cmake out-of-source builds
  configure)

    echo "Use recent CMake v3.26.3"
    export PATH=$HOME/opt/cmake/3.26.3/bin:$PATH
    # Make current environment variables available to subsequent steps, ctest
    echo "PATH=$PATH" >> $GITHUB_ENV

    QMC_DATA_DIR=/scratch/ci/QMC_DATA_FULL

    # Using 1.74 to avoid the > 1.75 error: use of undeclared identifier 'noinline'; did you mean 'inline'?
    # caused by LLVM + GCC libstdc++ mismatch
    BOOST_DIR=$HOME/opt/spack/linux-rhel9-cascadelake/gcc-9.4.0/boost-1.74.0-gdhlc5uynyw5un6mniss7nfjdyqqjd7p

    if [ -d ${GITHUB_WORKSPACE}/../qmcpack-build-2 ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build-2, removing"
      rm -fr ${GITHUB_WORKSPACE}/../qmcpack-build-2
    fi

    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build-2"
    cd ${GITHUB_WORKSPACE}/.. && mkdir qmcpack-build-2 && cd qmcpack-build-2

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
        echo 'Configure for mixed precision build -DQMC_MIXED_PRECISION=0'
        IS_MIXED_PRECISION=0
      ;;
    esac

    case "${GH_JOBNAME}" in

      *"V100-Clang16-MPI-CUDA-AFQMC-Offload"*)
        echo "Configure for building with CUDA and AFQMC using LLVM OpenMP offload"

        LLVM_DIR=$HOME/opt/llvm/18.1.2

        echo "Set PATHs to CUDA"
        export PATH=/usr/local/cuda/bin:$PATH

        export OMPI_CC=$LLVM_DIR/bin/clang
        export OMPI_CXX=$LLVM_DIR/bin/clang++

        # Make current environment variables available to subsequent steps
        echo "PATH=$PATH" >> $GITHUB_ENV
        echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $GITHUB_ENV
        echo "OMPI_CC=$OMPI_CC" >> $GITHUB_ENV
        echo "OMPI_CXX=$OMPI_CXX" >> $GITHUB_ENV

        # Confirm that CUDA and GCC toolchain get picked up by the compiler as expected
        $OMPI_CXX -v

        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBOOST_ROOT=$BOOST_DIR \
              -DBUILD_AFQMC=ON \
              -DQMC_GPU="cuda;openmp" \
              -DQMC_GPU_ARCHS=sm_70 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE} || cmake_exit_code=$?
      ;;

      *"V100-GCC11-MPI-CUDA"*)
        echo "Configure for building with CUDA and" \
             "GCC11 system compiler"

        echo "Set PATHs to CUDA"
        export PATH=/usr/local/cuda/bin:$PATH

        # Make current environment variables available to subsequent steps
        echo "PATH=$PATH" >> $GITHUB_ENV
        echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $GITHUB_ENV

        cmake -GNinja \
              -DCMAKE_C_COMPILER=/usr/lib64/openmpi/bin/mpicc \
              -DCMAKE_CXX_COMPILER=/usr/lib64/openmpi/bin/mpicxx \
              -DMPIEXEC_EXECUTABLE=/usr/lib64/openmpi/bin/mpirun \
              -DBOOST_ROOT=$BOOST_DIR \
              -DQMC_GPU=cuda \
              -DQMC_GPU_ARCHS=sm_70 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE} || cmake_exit_code=$?
      ;;
    esac
    if [ $cmake_exit_code -ne 0 ]; then
      # for debugging purpose
      if [ -f CMakeFiles/CMakeError.log ]; then cat CMakeFiles/CMakeError.log; fi
      if [ -f CMakeFiles/CMakeOutput.log ]; then cat CMakeFiles/CMakeOutput.log; fi
    fi
    exit $cmake_exit_code
    ;;

  build)
    # Verify nvcc
    which nvcc
    cd ${GITHUB_WORKSPACE}/../qmcpack-build-2
    ninja
    ;;

  test)
    echo "Enabling OpenMPI oversubscription"
    export OMPI_MCA_rmaps_base_oversubscribe=1
    export OMPI_MCA_hwloc_base_binding_policy=none
    echo "Set the management layer to ucx"
    export OMPI_MCA_pml=ucx
    # Avoid polluting the stderr output with libfabric error message
    export OMPI_MCA_btl=self
    # Clang helper threads used by target nowait is very broken. Disable this feature
    export LIBOMP_USE_HIDDEN_HELPER_TASK=0
    export OMP_TARGET_OFFLOAD=mandatory

    echo "Running deterministic tests"
    cd ${GITHUB_WORKSPACE}/../qmcpack-build-2
    if [ -f ./bin/qmcpack ]; then ldd ./bin/qmcpack; fi
    if [ -f ./bin/qmcpack_complex ]; then ldd ./bin/qmcpack_complex; fi
    ctest --output-on-failure -E ppconvert -L deterministic -j 32
    ;;

esac
