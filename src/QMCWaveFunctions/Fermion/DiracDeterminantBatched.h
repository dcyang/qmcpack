//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminantBatched.h
 * @brief Declaration of DiracDeterminantBatched with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H
#define QMCPLUSPLUS_DIRACDETERMINANTBATCHED_H

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "WaveFunctionTypes.hpp"
#include "type_traits/complex_help.hpp"
#include "QMCWaveFunctions/Fermion/DelayedUpdateBatched.h"
#include "DiracMatrixInverter.hpp"

namespace qmcplusplus
{

//forward declaration
class TWFFastDerivWrapper;

template<PlatformKind PL, typename VT>
struct UpdateEngineSelector;

template<typename VT>
struct UpdateEngineSelector<PlatformKind::OMPTARGET, VT>
{
  using Engine = DelayedUpdateBatched<PlatformKind::OMPTARGET, VT>;
};

#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
template<typename VT>
struct UpdateEngineSelector<PlatformKind::CUDA, VT>
{
  using Engine = DelayedUpdateBatched<PlatformKind::CUDA, VT>;
};
#endif

#if defined(ENABLE_SYCL) && defined(ENABLE_OFFLOAD)
template<typename VT>
struct UpdateEngineSelector<PlatformKind::SYCL, VT>
{
  using Engine = DelayedUpdateBatched<PlatformKind::SYCL, VT>;
};
#endif

template<PlatformKind PL, typename VT, typename FPVT>
class DiracDeterminantBatched : public DiracDeterminantBase
{
public:
  using UpdateEngine  = typename UpdateEngineSelector<PL, VT>::Engine;
  using WFT           = WaveFunctionTypes<VT, FPVT>;
  using Value         = typename WFT::Value;
  using FullPrecValue = typename WFT::FullPrecValue;
  using PsiValue      = typename WFT::PsiValue;
  using LogValue      = typename WFT::LogValue;
  using Grad          = typename WFT::Grad;
  using Hess          = typename WFT::Hess;
  using Real          = typename WFT::Real;
  using FullPrecGrad  = TinyVector<FullPrecValue, DIM>;

  // the understanding of dual memory space needs to follow UpdateEngine
  template<typename DT>
  using DualVector = Vector<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using DualMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  using DualVGLVector = VectorSoaContainer<Value, DIM + 2, OffloadPinnedAllocator<Value>>;

  using OffloadMWVGLArray = typename SPOSet::OffloadMWVGLArray;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   *@param last index of last particle
   *@param ndelay delayed update rank
   */
  DiracDeterminantBatched(SPOSet& phi,
                          int first,
                          int last,
                          int ndelay                          = 1,
                          DetMatInvertor matrix_inverter_kind = DetMatInvertor::ACCEL);

  // copy constructor and assign operator disabled
  DiracDeterminantBatched(const DiracDeterminantBatched& s)            = delete;
  DiracDeterminantBatched& operator=(const DiracDeterminantBatched& s) = delete;

  std::string getClassName() const override { return "DiracDeterminant"; }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           Vector<Value>& dlogpsi,
                           Vector<Value>& dhpsioverpsi) override;

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValue ratio(ParticleSet& P, int iat) override;

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios) const override;

  /** compute multiple ratios for a particle move
   */
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<Value>& ratios) override;

  void evaluateSpinorRatios(const VirtualParticleSet& VP,
                            const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                            std::vector<Value>& ratios) override;

  void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                         std::vector<std::vector<Value>>& ratios) const override;

  void mw_evaluateSpinorRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                               const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                               const RefVector<std::pair<ValueVector, ValueVector>>& spinor_multiplier_list,
                               std::vector<std::vector<Value>>& ratios) const override;

  void evaluateDerivRatios(const VirtualParticleSet& VP,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& ratios,
                           Matrix<ValueType>& dratios) override;

  void evaluateSpinorDerivRatios(const VirtualParticleSet& VP,
                                 const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                 const opt_variables_type& optvars,
                                 std::vector<ValueType>& ratios,
                                 Matrix<ValueType>& dratios) override;

  PsiValue ratioGrad(ParticleSet& P, int iat, Grad& grad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<Grad>& grad_new) const override;

  void mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValue>& ratios,
                            std::vector<Grad>& grad_new,
                            std::vector<ComplexType>& spingrad_new) const override;

  PsiValue ratioGradWithSpin(ParticleSet& P, int iat, Grad& grad_iat, ComplexType& spingrad) override;

  Grad evalGrad(ParticleSet& P, int iat) override;

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<Grad>& grad_now) const override;

  Grad evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override;

  void mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<Grad>& grad_now,
                           std::vector<ComplexType>& spingrad_now) const override;

  /** \todo would be great to have docs.
   *  Note: Can result in substantial CPU memory allocation on first call.
   *  31 * n^2 * sizeof(Value) bytes per DDB
   */
  Grad evalGradSource(ParticleSet& P, ParticleSet& source, int iat) override;

  Grad evalGradSource(ParticleSet& P,
                      ParticleSet& source,
                      int iat,
                      TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                      TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override;

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay = false) const override;

  /** complete any left over determinant matrix updates.
   * Usually this is the end of pbyp moves for a given spin of electrons
   * The psiM, dpsiM, d2psiM should be up-to-date on both device and host sides.
   */
  void completeUpdates() override;

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  /** evaluate from scratch pretty much everything for this single walker determinant
   *
   *  return of the log of the dirac determinant is the least of what it does.
   *
   *  call to generate valid initial state for determinant and when you
   *  suspect psiMinv or other state variables may have picked up error.
   */
  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian>& L_list) const override;

  void recompute(const ParticleSet& P) override;

  /** Does a phi_.mw_evaluate_notranspose then mw_invertPsiM over a set of
   *  elements filtered based on the recompute mask.
   *
   */
  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;

  LogValue evaluateGL(const ParticleSet& P,
                      ParticleSet::ParticleGradient& G,
                      ParticleSet::ParticleLaplacian& L,
                      bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                     bool fromscratch) const override;

  void evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi) override;

  void createResource(ResourceCollection& collection) const override;
  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;
  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  void registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const override;
  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  std::unique_ptr<DiracDeterminantBase> makeCopy(SPOSet& phi) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<Value>& ratios) override;

  const auto& get_psiMinv() const { return psiMinv_; }

private:
  /** @defgroup LegacySingleData
   *  @brief    Single Walker Data Members of Legacy OO design
   *            High and flexible throughput of walkers requires would ideally separate
   *            walker data which should be "SoA" and functions over it i.e. leave behind
   *            the OO pattern of a single set of data and functions on it.
   *  
   *  @ingroup LegacySingleData
   *  @{
   */
  /* inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
   * Only NumOrbitals x NumOrbitals subblock has meaningful data
   * The number of rows is equal to NumOrbitals
   * The number of columns in each row is padded to a multiple of QMC_SIMD_ALIGNMENT
   */
  DualMatrix<Value> psiMinv_;
  /// fused memory for psiM, dpsiM and d2psiM. [5][norb*norb]
  DualVGLVector psiM_vgl;
  /** psiM(j,i) \f$= \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
   *  Only this one is Dual since only psiM is a Dual argument
   *  in the det_engine single walker API.
   */
  DualMatrix<Value> psiM_temp;
  Matrix<Value> psiM_host;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  Matrix<Grad> dpsiM;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$. partial memory view of psiM_vgl
  Matrix<Value> d2psiM;

  /// Used for force computations
  Matrix<Grad> grad_source_psiM, grad_lapl_source_psiM;
  Matrix<Hess> grad_grad_source_psiM;

  Matrix<Grad> phi_alpha_Minv, grad_phi_Minv;
  Matrix<Value> lapl_phi_Minv;
  Matrix<Hess> grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  DualVector<Value> psiV;
  Vector<Value> psiV_host_view;
  DualVector<Grad> dpsiV;
  Vector<Grad> dpsiV_host_view;
  DualVector<Value> d2psiV;
  Vector<Value> d2psiV_host_view;
  DualVector<Value> dspin_psiV;
  Vector<Value> dspin_psiV_host_view;

  /// psi(r')/psi(r) during a PbyP move
  PsiValue curRatio;
  /**@}*/

  struct DiracDeterminantBatchedMultiWalkerResource;
  ResourceHandle<DiracDeterminantBatchedMultiWalkerResource> mw_res_handle_;

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  /// Delayed update engine 1 per walker.
  UpdateEngine det_engine_;

  /// slow but doesn't consume device memory
  DiracMatrix<FullPrecValue> host_inverter_;

  /// matrix inversion engine this a crowd scope resource and only the leader engine gets it
  ResourceHandle<DiracMatrixInverter<FPVT, VT>> accel_inverter_;

  /// compute G and L assuming psiMinv, dpsiM, d2psiM are ready for use
  void computeGL(ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L) const;

  /// single invert logdetT(psiM)
  /// as a side effect this->log_value_ gets the log determinant of logdetT
  void invertPsiM(const DualMatrix<Value>& psiM, DualMatrix<Value>& psiMinv);

  /** Inverts and finds log det for a batch of matrices
   *
   *  Right now this takes filtered lists and the full recompute mask. It passes
   *  all these elements down to the det_engine_.
   *  This allows minimal change for implementation code while establishing the API
   *  I'd prefer for direct inversion.
   *  I think the det_engine_ mw method  should receive the complete
   *  list of walker elements and the implementation should decide what to do re
   *  the compute mask. See future PR for those changes, or drop of compute_mask argument.
   */
  static void mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVector<const DualMatrix<Value>>& logdetT_list,
                            const RefVector<DualMatrix<Value>>& a_inv_lis);

  /// Resize all temporary arrays required for force computation.
  void resizeScratchObjectsForIonDerivs();

  /// maximal number of delayed updates
  const int ndelay_;

  /// selected scheme for inversion with walker batching
  const DetMatInvertor matrix_inverter_kind_;

  /// timers
  NewTimer &D2HTimer, &H2DTimer;
};

extern template class DiracDeterminantBatched<PlatformKind::OMPTARGET,
                                              QMCTraits::ValueType,
                                              QMCTraits::QTFull::ValueType>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
extern template class DiracDeterminantBatched<PlatformKind::CUDA, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif
#if defined(ENABLE_SYCL) && defined(ENABLE_OFFLOAD)
extern template class DiracDeterminantBatched<PlatformKind::SYCL, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif


} // namespace qmcplusplus
#endif
