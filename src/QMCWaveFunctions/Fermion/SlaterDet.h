//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"

namespace qmcplusplus
{
class TWFFastDerivWrapper;

class SlaterDet : public WaveFunctionComponent
{
public:
  using Determinant_t = DiracDeterminantBase;

  /**  constructor
   * @param targetPtcl target Particleset
   */
  SlaterDet(ParticleSet& targetPtcl,
            std::vector<std::unique_ptr<SPOSet>>&& sposets,
            std::vector<std::unique_ptr<Determinant_t>>&& dets,
            const std::string& class_name = "SlaterDet");

  ///destructor
  ~SlaterDet() override;

  std::string getClassName() const override { return "SlaterDet"; }

  bool isFermionic() const final { return true; }
  bool isOptimizable() const override;

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override;

  void checkOutVariables(const opt_variables_type& active) override;

  void registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const override;

  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian>& L_list) const override;

  LogValue evaluateGL(const ParticleSet& P,
                      ParticleSet::ParticleGradient& G,
                      ParticleSet::ParticleLaplacian& L,
                      bool fromscratch) override;

  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                     bool fromscratch) const override;

  void recompute(const ParticleSet& P) override;

  void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    const std::vector<bool>& recompute) const override;

  void evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  void createResource(ResourceCollection& collection) const override;

  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;

  inline void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override
  {
    return Dets[getDetID(VP.refPtcl)]->evaluateRatios(VP, ratios);
  }

  inline void evaluateSpinorRatios(const VirtualParticleSet& VP,
                                   const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                   std::vector<ValueType>& ratios) override
  {
    return Dets[getDetID(VP.refPtcl)]->evaluateSpinorRatios(VP, spinor_multiplier, ratios);
  }

  void evaluateDerivRatios(const VirtualParticleSet& VP,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& ratios,
                           Matrix<ValueType>& dratios) override;

  void evaluateSpinorDerivRatios(const VirtualParticleSet& VP,
                                 const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                 const opt_variables_type& optvars,
                                 std::vector<ValueType>& ratios,
                                 Matrix<ValueType>& dratios) override;

  inline void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                std::vector<std::vector<ValueType>>& ratios) const override
  {
    if (wfc_list.size())
    {
      // assuming all the VP.refPtcl are identical
      const int det_id = getDetID(vp_list[0].refPtcl);
      Dets[det_id]->mw_evaluateRatios(extract_DetRef_list(wfc_list, det_id), vp_list, ratios);
    }
  }

  inline void mw_evaluateSpinorRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                      const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                      const RefVector<std::pair<ValueVector, ValueVector>>& spinor_multiplier_list,
                                      std::vector<std::vector<ValueType>>& ratios) const final
  {
    if (wfc_list.size())
    {
      // assuming all the VP.refPtcl are identical
      const int det_id = getDetID(vp_list[0].refPtcl);
      Dets[det_id]->mw_evaluateSpinorRatios(extract_DetRef_list(wfc_list, det_id), vp_list, spinor_multiplier_list,
                                            ratios);
    }
  }

  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  PsiValue ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat) override;

  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<GradType>& grad_now) const override;

  void mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValue>& ratios,
                            std::vector<GradType>& grad_now,
                            std::vector<ComplexType>& spingrad_now) const override;

  GradType evalGrad(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->evalGrad(P, iat); }

  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    return Dets[getDetID(iat)]->evalGradWithSpin(P, iat, spingrad);
  }

  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_evalGrad(extract_DetRef_list(wfc_list, det_id), p_list, iat, grad_now);
  }

  void mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<GradType>& grad_now,
                           std::vector<ComplexType>& spingrad_now) const override;

  GradType evalGradSource(ParticleSet& P, ParticleSet& src, int iat) override
  {
    GradType G = GradType();
    for (int iz = 0; iz < Dets.size(); iz++)
      G += Dets[iz]->evalGradSource(P, src, iat);
    return G;
  }

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& src,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override
  {
    GradType G = GradType();
    for (int iz = 0; iz < Dets.size(); iz++)
      G += Dets[iz]->evalGradSource(P, src, iat, grad_grad, lapl_grad);
    return G;
  }

  inline void restore(int iat) override { return Dets[getDetID(iat)]->restore(iat); }

  inline void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    Dets[getDetID(iat)]->acceptMove(P, iat, safe_to_delay);

    log_value_ = 0.0;
    for (int i = 0; i < Dets.size(); ++i)
      log_value_ += Dets[i]->get_log_value();
  }

  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay = false) const override
  {
    constexpr LogValue czero(0);

    // This log_value_ is in the slater determinant, it's still around but not consistent anymore with the
    // sum of the log_values in its determinants.  Caching the state seems like a bad call, but the wfc base class
    // having log_value_ as a data member asks for this sort of consistency issue when wfc can contain wfc.
    for (int iw = 0; iw < wfc_list.size(); iw++)
      if (isAccepted[iw])
        wfc_list.getCastedElement<SlaterDet>(iw).log_value_ = czero;

    for (int i = 0; i < Dets.size(); ++i)
    {
      const auto Det_list(extract_DetRef_list(wfc_list, i));

      if (i == getDetID(iat))
        Dets[i]->mw_accept_rejectMove(Det_list, p_list, iat, isAccepted, safe_to_delay);

      for (int iw = 0; iw < wfc_list.size(); iw++)
        if (isAccepted[iw])
          wfc_list.getCastedElement<SlaterDet>(iw).log_value_ += Det_list[iw].get_log_value();
    }
  }

  void completeUpdates() override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->completeUpdates();
  }

  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override
  {
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->mw_completeUpdates(extract_DetRef_list(wfc_list, i));
  }

  inline PsiValue ratio(ParticleSet& P, int iat) override { return Dets[getDetID(iat)]->ratio(P, iat); }

  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios) const override
  {
    const int det_id = getDetID(iat);
    Dets[det_id]->mw_calcRatio(extract_DetRef_list(wfc_list, det_id), p_list, iat, ratios);
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override
  {
    // Now add on contribution from each determinant to the derivatives
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
  }

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, Vector<ValueType>& dlogpsi) override
  {
    // Now add on contribution from each determinant to the derivatives
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->evaluateDerivativesWF(P, active, dlogpsi);
  }

  ///return the total number of Dirac determinants
  inline int getNumDets() const { return Dets.size(); }
  ///return the i-th determinant
  inline auto& getDet(const int i) { return *Dets[i]; }
  ///return the sposet of the i-th determinant
  SPOSet& getPhi(int i = 0) { return Dets[i]->getPhi(); }

private:
  //get Det ID
  inline int getDetID(const int iat) const
  {
    int id = 0;
    while (iat > Last[id])
      id++;
    return id;
  }

  // helper function for extracting a list of WaveFunctionComponent from a list of TrialWaveFunction
  RefVectorWithLeader<WaveFunctionComponent> extract_DetRef_list(
      const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
      int det_id) const
  {
    RefVectorWithLeader<WaveFunctionComponent> Det_list(*wfc_list.getCastedLeader<SlaterDet>().Dets[det_id]);
    Det_list.reserve(wfc_list.size());
    for (WaveFunctionComponent& wfc : wfc_list)
      Det_list.push_back(*static_cast<SlaterDet&>(wfc).Dets[det_id]);
    return Det_list;
  }

  ///the last particle of each group
  std::vector<int> Last;

  ///container for the unique SPOSets
  const std::vector<std::unique_ptr<SPOSet>> sposets_;

  ///container for the DiracDeterminants
  const std::vector<std::unique_ptr<Determinant_t>> Dets;
};
} // namespace qmcplusplus
#endif
