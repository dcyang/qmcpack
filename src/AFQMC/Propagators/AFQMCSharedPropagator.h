//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SHAREDPROPAGATOR_H
#define QMCPLUSPLUS_AFQMC_SHAREDPROPAGATOR_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"

namespace qmcplusplus
{

namespace afqmc
{

/*
 * AFQMC propagator using shared memory. Data structures are in shared memory, 
 * but not distributed over multiple nodes. 
 */
class AFQMCSharedPropagator: public AFQMCInfo
{
  protected:

  using CVector = boost::multi::array<ComplexType,1>;  
  using CVector_ref = boost::multi::array_ref<ComplexType,1>;  
  using CMatrix = boost::multi::array<ComplexType,2>;  
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2>;  
  using shmCVector = ComplexVector<shared_allocator<ComplexType>>; 

  public:

    AFQMCSharedPropagator(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, 
                          Wavefunction& wfn_, CMatrix&& h1_, CVector&& vmf_, 
                          RandomGenerator_t* r):
            AFQMCInfo(info),TG(tg_),wfn(wfn_),
            H1(std::move(h1_)),
            P1(P1Type(tp_ul_ul{0,0},tp_ul_ul{0,0},0,boost::mpi3::intranode::allocator<ComplexType>(tg_.TG_local()))),
            vMF(std::move(vmf_)),
            rng(r),
            SDetOp(2*NMO,NAEA+NAEB), // safe for now, since I don't know walker_type
            TSM({2*NMO,NAEA+NAEB}), // safe for now, since I don't know walker_type
            shmbuff(extensions<1u>{1},shared_allocator<ComplexType>{TG.TG_local()}),
            local_group_comm(),
            last_nextra(-1),
            last_task_index(-1),
            old_dt(-123456.789),
            order(6)
    {
      transposed_vHS_ = wfn.transposed_vHS();
      transposed_G_ = wfn.transposed_G_for_vbias();
      if(not transposed_vHS_) local_vHS.reextent({NMO,NMO});
      parse(cur);  
    }

    ~AFQMCSharedPropagator() {}

    AFQMCSharedPropagator(AFQMCSharedPropagator const& other) = delete;
    AFQMCSharedPropagator& operator=(AFQMCSharedPropagator const& other) = delete;
    AFQMCSharedPropagator(AFQMCSharedPropagator&& other) = default;
    AFQMCSharedPropagator& operator=(AFQMCSharedPropagator&& other) = default;

    template<class WlkSet>
    void Propagate(int steps, WlkSet& wset, RealType E1,
                   RealType dt, int fix_bias=1) { 
      int nblk = steps/fix_bias;
      int nextra = steps%fix_bias;
      for(int i=0; i<nblk; i++) {
        step(fix_bias,wset,E1,dt);
      }
      if(nextra>0) 
        step(nextra,wset,E1,dt);
      TG.local_barrier();
    }

    // reset shared memory buffers
    // useful when the current buffers use too much memory (e.g. reducing steps in future calls)
    void reset() { shmbuff.reextent(extensions<1u>{0}); }

    bool hybrid_propagation() { return hybrid; }

    bool free_propagation() { return free_projection; }

  protected: 

    TaskGroup_& TG;

    Wavefunction& wfn;

    // P1 = exp(-0.5*dt*H1), so H1 includes terms from MF substraction 
    //                       and the exchange term from the cholesky decomposition (e.g. vn0) 
    CMatrix H1;

    P1Type P1;

    RandomGenerator_t* rng;

    SlaterDetOperations<ComplexType> SDetOp;

    shmCVector shmbuff;    

    shared_communicator local_group_comm;

    RealType old_dt;
    int last_nextra;
    int last_task_index;
    int order;

    RealType vbias_bound;

    // type of propagation
    bool free_projection;
    bool hybrid;
    bool importance_sampling;
    bool apply_constrain;

    bool transposed_vHS_;
    bool transposed_G_;

    // used in propagator step 
    CMatrix local_vHS;  

    CVector new_overlaps;  
    CMatrix new_energies;  

    CMatrix MFfactor;  
    CMatrix hybrid_weight;  

    CVector vMF;  
    // Temporary for propagating with constructed B matrix.
    CMatrix TSM;
 
    template<class WlkSet>
    void step(int steps, WlkSet& wset, RealType E1, RealType dt);

    void assemble_X(size_t nsteps, size_t nwalk, RealType sqrtdt,
                    CMatrix_ref& X, CMatrix_ref & vbias,
                    CMatrix_ref& MF, CMatrix_ref& HW, bool addRAND=true);

    void reset_nextra(int nextra); 

    void parse(xmlNodePtr cur);

    template<class WSet>
    void apply_propagators(WSet& wset, int ni, int tk0, int tkN, int ntask_total_serial,
                           boost::multi::array_ref<ComplexType,3>& vHS3D);

    template<class WSet>
    void apply_propagators_construct_propagator(WSet& wset, int ni, int tk0, int tkN, int ntask_total_serial,
                                                boost::multi::array_ref<ComplexType,3>& vHS3D);

    ComplexType apply_bound_vbias(ComplexType v, RealType sqrtdt)
    {
      return (std::abs(v)>vbias_bound*sqrtdt)?
                (v/(std::abs(v)/static_cast<ValueType>(vbias_bound*sqrtdt))):(v);
    }
  
};

}

}

#include "AFQMCSharedPropagator.icc"

#endif

