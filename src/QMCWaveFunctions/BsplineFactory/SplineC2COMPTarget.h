//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file SplineC2COMPTarget.h
 *
 * class to handle complex splines to complex orbitals with splines of arbitrary precision
 * splines storage and computation is offloaded to accelerators using OpenMP target
 */
#ifndef QMCPLUSPLUS_SPLINE_C2C_OMPTARGET_H
#define QMCPLUSPLUS_SPLINE_C2C_OMPTARGET_H

#include <memory>
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBsplineOffload.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "Utilities/FairDivide.h"
#include "Utilities/TimerManager.h"
#include <ResourceHandle.h>
#include "SplineOMPTargetMultiWalkerMem.h"

namespace qmcplusplus
{
/** class to match std::complex<ST> spline with BsplineSet::ValueType (complex) SPOs with OpenMP offload
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of phase vectors
 * The internal storage of complex spline coefficients uses double sized real arrays of ST type, aligned and padded.
 * All the output orbitals are complex.
 */
template<typename ST>
class SplineC2COMPTarget : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using ComplexT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector;
  using BsplineSet::GradVector;
  using BsplineSet::HessVector;
  using BsplineSet::ValueVector;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

  template<typename DT>
  using OffloadVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using OffloadPosVector = VectorSoaContainer<DT, 3, OffloadAllocator<DT>>;

private:
  /// timer for offload portion
  NewTimer& offload_timer_;
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///multi bspline set
  std::shared_ptr<MultiBsplineBase<ST>> SplineInst;

  ///Copy of original splines for orbital rotation. Only need these on host
  std::shared_ptr<std::vector<ST>> coef_copy_;

  std::shared_ptr<OffloadVector<ST>> mKK;
  std::shared_ptr<OffloadPosVector<ST>> myKcart;
  std::shared_ptr<OffloadVector<ST>> GGt_offload;
  std::shared_ptr<OffloadVector<ST>> PrimLattice_G_offload;

  ResourceHandle<SplineOMPTargetMultiWalkerMem<ST, ComplexT>> mw_mem_handle_;

  ///team private ratios for reduction, numVP x numTeams
  Matrix<ComplexT, OffloadPinnedAllocator<ComplexT>> ratios_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<ComplexT, OffloadPinnedAllocator<ComplexT>> results_scratch;
  ///psiinv and position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ComplexT, OffloadPinnedAllocator<ComplexT>> psiinv_pos_copy;
  ///position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ST, OffloadPinnedAllocator<ST>> multi_pos_copy;

  void evaluateVGLMultiPos(const Vector<ST, OffloadPinnedAllocator<ST>>& multi_pos_copy,
                           Vector<ST, OffloadPinnedAllocator<ST>>& offload_scratch,
                           Vector<ComplexT, OffloadPinnedAllocator<ComplexT>>& results_scratch,
                           const RefVector<ValueVector>& psi_v_list,
                           const RefVector<GradVector>& dpsi_v_list,
                           const RefVector<ValueVector>& d2psi_v_list) const;

protected:
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2COMPTarget(const std::string& my_name, bool use_offload = true)
      : BsplineSet(my_name),
        offload_timer_(createGlobalTimer("SplineC2COMPTarget::offload", timer_level_fine)),
        GGt_offload(std::make_shared<OffloadVector<ST>>(9)),
        PrimLattice_G_offload(std::make_shared<OffloadVector<ST>>(9))
  {}

  SplineC2COMPTarget(const SplineC2COMPTarget& in);

  virtual std::string getClassName() const override { return "SplineC2COMPTarget"; }
  virtual std::string getKeyword() const override { return "SplineC2C"; }
  bool isComplex() const override { return true; };
  virtual bool isOMPoffload() const override { return true; }

  void createResource(ResourceCollection& collection) const override
  {
    auto resource_index = collection.addResource(std::make_unique<SplineOMPTargetMultiWalkerMem<ST, ComplexT>>());
  }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override
  {
    assert(this == &spo_list.getLeader());
    auto& phi_leader          = spo_list.getCastedLeader<SplineC2COMPTarget<ST>>();
    phi_leader.mw_mem_handle_ = collection.lendResource<SplineOMPTargetMultiWalkerMem<ST, ComplexT>>();
  }

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override
  {
    assert(this == &spo_list.getLeader());
    auto& phi_leader = spo_list.getCastedLeader<SplineC2COMPTarget<ST>>();
    collection.takebackResource(phi_leader.mw_mem_handle_);
  }

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<SplineC2COMPTarget>(*this); }

  bool isRotationSupported() const override { return true; }

  ///store copy of spline coefficients for obrital rotation
  void storeParamsBeforeRotation() override;

  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) override;

  inline void resizeStorage(size_t n) override
  {
    init_base(n);
    size_t npad = getAlignedSize<ST>(2 * n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
    mygH.resize(npad);
  }

  void bcast_tables(Communicate* comm) { chunked_bcast(comm, SplineInst->getSplinePtr()); }

  void gather_tables(Communicate* comm)
  {
    if (comm->size() == 1)
      return;
    const int Nbands      = kPoints.size();
    const int Nbandgroups = comm->size();
    offset.resize(Nbandgroups + 1, 0);
    FairDivideLow(Nbands, Nbandgroups, offset);

    for (size_t ib = 0; ib < offset.size(); ib++)
      offset[ib] *= 2;
    gatherv(comm, SplineInst->getSplinePtr(), SplineInst->getSplinePtr()->z_stride, offset);
  }

  template<typename BCT>
  void create_spline(const Ugrid xyz_g[3], const BCT& xyz_bc)
  {
    resize_kpoints();
    SplineInst = std::make_shared<MultiBsplineOffload<ST>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  /// this routine can not be called from threaded region
  void finalizeConstruction() override
  {
    SplineInst->finalize();
    // transfer static data to GPU
    mKK->updateTo();
    myKcart->updateTo();
    for (uint32_t i = 0; i < 9; i++)
    {
      (*GGt_offload)[i]           = GGt[i];
      (*PrimLattice_G_offload)[i] = PrimLattice.G[i];
    }
    PrimLattice_G_offload->updateTo();
    GGt_offload->updateTo();
  }

  inline void flush_zero() { SplineInst->flush_zero(); }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
    const size_t nk = kPoints.size();
    mKK             = std::make_shared<OffloadVector<ST>>(nk);
    myKcart         = std::make_shared<OffloadPosVector<ST>>(nk);
    for (size_t i = 0; i < nk; ++i)
    {
      (*mKK)[i]     = -dot(kPoints[i], kPoints[i]);
      (*myKcart)(i) = kPoints[i];
    }
  }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level);

  bool read_splines(hdf_archive& h5f);

  bool write_splines(hdf_archive& h5f);

  void assign_v(const PointType& r, const vContainer_type& myV, ValueVector& psi, int first, int last) const;

  virtual void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override;

  virtual void evaluateDetRatios(const VirtualParticleSet& VP,
                                 ValueVector& psi,
                                 const ValueVector& psiinv,
                                 std::vector<ValueType>& ratios) override;

  virtual void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                    const RefVector<ValueVector>& psi_list,
                                    const std::vector<const ValueType*>& invRow_ptr_list,
                                    std::vector<std::vector<ValueType>>& ratios_list) const override;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  virtual void evaluateVGL(const ParticleSet& P,
                           const int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           ValueVector& d2psi) override;

  virtual void mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& sa_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector>& psi_v_list,
                              const RefVector<GradVector>& dpsi_v_list,
                              const RefVector<ValueVector>& d2psi_v_list) const override;

  virtual void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                              const RefVectorWithLeader<ParticleSet>& P_list,
                                              int iat,
                                              const std::vector<const ValueType*>& invRow_ptr_list,
                                              OffloadMWVGLArray& phi_vgl_v,
                                              std::vector<ValueType>& ratios,
                                              std::vector<GradType>& grads) const override;

  void assign_vgh(const PointType& r,
                  ValueVector& psi,
                  GradVector& dpsi,
                  HessVector& grad_grad_psi,
                  int first,
                  int last) const;

  virtual void evaluateVGH(const ParticleSet& P,
                           const int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& grad_grad_psi) override;

  void assign_vghgh(const PointType& r,
                    ValueVector& psi,
                    GradVector& dpsi,
                    HessVector& grad_grad_psi,
                    GGGVector& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  virtual void evaluateVGHGH(const ParticleSet& P,
                             const int iat,
                             ValueVector& psi,
                             GradVector& dpsi,
                             HessVector& grad_grad_psi,
                             GGGVector& grad_grad_grad_psi) override;

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix& logdet,
                                    GradMatrix& dlogdet,
                                    ValueMatrix& d2logdet) override;

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend struct BsplineReader;
};

extern template class SplineC2COMPTarget<float>;
extern template class SplineC2COMPTarget<double>;

} // namespace qmcplusplus
#endif
