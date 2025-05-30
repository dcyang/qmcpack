.. _chap:features:

Features of QMCPACK
===================

Note that besides direct use, most features are also available via Nexus, an advanced workflow tool to automate all aspects of QMC
calculation from initial DFT calculations through to final analysis. Use of Nexus is highly recommended for research calculations
due to the greater ease of use and increased reproducibility.

Real-space Monte Carlo
----------------------

The following list contains the main production-level features of QMCPACK for real-space Monte Carlo. If you do not see a specific
feature that you are interested in, check the remainder of this manual or ask if that specific feature can be made available.

-  Variational Monte Carlo (VMC).

-  Diffusion Monte Carlo (DMC).

-  Reptation Monte Carlo.

-  Single and multideterminant Slater Jastrow wavefunctions.

-  Wavefunction updates using optimized multideterminant algorithm of
   Clark et al.

-  Backflow wavefunctions.

-  One, two, and three-body Jastrow factors.

-  Excited state calculations via flexible occupancy assignment of
   Slater determinants.

-  All electron and nonlocal pseudopotential calculations.

-  Casula T-moves for variational evaluation of nonlocal
   pseudopotentials (non-size-consistent and size-consistent variants).

-  Spin-orbit coupling from relativistic pseudopotentials following the 
   approach of Melton, Bennett, and Mitas.

-  Support for twist boundary conditions and calculations on metals.

-  Wavefunction optimization using the “linear method” of Umrigar and
   coworkers, with an arbitrary mix of variance and energy in the objective
   function.

-  Blocked, low memory adaptive shift optimizer of Zhao and Neuscamman.

-  Gaussian, Slater, plane-wave, and real-space spline basis sets for
   orbitals.

-  Interface and conversion utilities for plane-wave wavefunctions from
   Quantum ESPRESSO (Plane-Wave Self-Consistent Field package [PWSCF]).

-  Interface and conversion utilities for Gaussian-basis wavefunctions
   from GAMESS, PySCF, and QP2. Many more are supported via the molden format and molden2qmc.

-  Easy extension and interfacing to other electronic structure codes
   via standardized XML and HDF5 inputs.

-  MPI parallelism, with scaling to millions of cores.

-  Fully threaded on CPUs using OpenMP.

-  Highly efficient vectorized CPU code tailored for modern architectures. :cite:`IPCC_SC17`

-  OpenMP-offload-based performance portable GPU implementation. Fully supports NVIDIA, AMD, and Intel GPUs.
   GPU and CPU execution can be mixed and matched.

-  Analysis tools for minimal environments (Perl only) through to
   Python-based environments with graphs produced via matplotlib (included with Nexus).

Auxiliary-Field Quantum Monte Carlo
-----------------------------------

The orbital-space Auxiliary-Field Quantum Monte Carlo (AFQMC) method is now also available in QMCPACK. The main input data are the
matrix elements of the Hamiltonian in a given single particle basis set, which must be produced from mean-field calculations such
as Hartree-Fock or density functional theory. A partial list of the current capabilities of the code follows. For a detailed
description of the available features, see  :ref:`afqmc`.

-  Phaseless AFQMC algorithm of Zhang et al. :cite:`PhysRevLett.90.136401`.

-  Very efficient GPU implementation for most features. 

-  “Hybrid" and “local energy" propagation schemes.

-  Hamiltonian matrix elements from (1) Molpro’s FCIDUMP format (which
   can be produced by Molpro, PySCF, and VASP) and (2) internal HDF5
   format produced by PySCF (see AFQMC section below).

-  AFQMC calculations with RHF (closed-shell doubly occupied), ROHF
   (open-shell doubly occupied), and UHF (spin polarized broken
   symmetry) symmetry.

-  Single and multideterminant trial wavefunctions. Multideterminant
   expansions with either orthogonal or nonorthogonal determinants.

-  Fast update scheme for orthogonal multideterminant expansions.

-  Distributed propagation algorithms for large systems. Enables
   calculations where data structures do not fit on a single node.

-  Complex implementation for PBC calculations with complex integrals.

-  Sparse representation of large matrices for reduced memory usage.

-  Mixed and back-propagated estimators.

-  Specialized implementation for solids with k-point symmetry (e.g.
   primitive unit cells with k-points).
