//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SlaterDet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

// for return types
using PsiValue = WaveFunctionComponent::PsiValue;

SlaterDet::SlaterDet(ParticleSet& targetPtcl,
                     std::vector<std::unique_ptr<SPOSet>>&& sposets,
                     std::vector<std::unique_ptr<Determinant_t>>&& dets,
                     const std::string& class_name)
    : sposets_(std::move(sposets)), Dets(std::move(dets))
{
  assert(Dets.size() == targetPtcl.groups());

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i) - 1;
}

///destructor
SlaterDet::~SlaterDet() = default;

bool SlaterDet::isOptimizable() const
{
  return std::any_of(sposets_.begin(), sposets_.end(), [](const auto& phi) { return phi->isOptimizable(); });
}

void SlaterDet::extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs)
{
  for (const auto& sposet : sposets_)
    sposet->extractOptimizableObjectRefs(opt_obj_refs);
}

void SlaterDet::checkOutVariables(const opt_variables_type& active)
{
  for (const auto& sposet : sposets_)
    if (sposet->isOptimizable())
      sposet->checkOutVariables(active);
}

PsiValue SlaterDet::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  return Dets[getDetID(iat)]->ratioGrad(P, iat, grad_iat);
}

PsiValue SlaterDet::ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat)
{
  return Dets[getDetID(iat)]->ratioGradWithSpin(P, iat, grad_iat, spingrad_iat);
}

void SlaterDet::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             int iat,
                             std::vector<PsiValue>& ratios,
                             std::vector<GradType>& grad_now) const
{
  const int det_id = getDetID(iat);
  Dets[det_id]->mw_ratioGrad(extract_DetRef_list(wfc_list, det_id), p_list, iat, ratios, grad_now);
}

void SlaterDet::mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                     const RefVectorWithLeader<ParticleSet>& p_list,
                                     int iat,
                                     std::vector<PsiValue>& ratios,
                                     std::vector<GradType>& grad_now,
                                     std::vector<ComplexType>& spingrad_now) const
{
  const int det_id = getDetID(iat);
  Dets[det_id]->mw_ratioGradWithSpin(extract_DetRef_list(wfc_list, det_id), p_list, iat, ratios, grad_now,
                                     spingrad_now);
}


void SlaterDet::mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    int iat,
                                    std::vector<GradType>& grad_now,
                                    std::vector<ComplexType>& spingrad_now) const
{
  const int det_id = getDetID(iat);
  Dets[det_id]->mw_evalGradWithSpin(extract_DetRef_list(wfc_list, det_id), p_list, iat, grad_now, spingrad_now);
}

void SlaterDet::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->evaluateRatiosAlltoOne(P, ratios);
}

void SlaterDet::evaluateDerivRatios(const VirtualParticleSet& VP,
                                    const opt_variables_type& optvars,
                                    std::vector<ValueType>& ratios,
                                    Matrix<ValueType>& dratios)
{
  return Dets[getDetID(VP.refPtcl)]->evaluateDerivRatios(VP, optvars, ratios, dratios);
}

void SlaterDet::evaluateSpinorDerivRatios(const VirtualParticleSet& VP,
                                          const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                          const opt_variables_type& optvars,
                                          std::vector<ValueType>& ratios,
                                          Matrix<ValueType>& dratios)
{
  return Dets[getDetID(VP.refPtcl)]->evaluateSpinorDerivRatios(VP, spinor_multiplier, optvars, ratios, dratios);
}

SlaterDet::LogValue SlaterDet::evaluateLog(const ParticleSet& P,
                                           ParticleSet::ParticleGradient& G,
                                           ParticleSet::ParticleLaplacian& L)
{
  log_value_ = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    log_value_ += Dets[i]->evaluateLog(P, G, L);
  return log_value_;
}

void SlaterDet::mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                               const RefVectorWithLeader<ParticleSet>& p_list,
                               const RefVector<ParticleSet::ParticleGradient>& G_list,
                               const RefVector<ParticleSet::ParticleLaplacian>& L_list) const
{
  constexpr LogValue czero(0);

  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list.getCastedElement<SlaterDet>(iw).log_value_ = czero;

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto Det_list(extract_DetRef_list(wfc_list, i));
    Dets[i]->mw_evaluateLog(Det_list, p_list, G_list, L_list);
    for (int iw = 0; iw < wfc_list.size(); iw++)
      wfc_list.getCastedElement<SlaterDet>(iw).log_value_ += Det_list[iw].get_log_value();
  }
}

SlaterDet::LogValue SlaterDet::evaluateGL(const ParticleSet& P,
                                          ParticleSet::ParticleGradient& G,
                                          ParticleSet::ParticleLaplacian& L,
                                          bool from_scratch)
{
  log_value_ = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    log_value_ += Dets[i]->evaluateGL(P, G, L, from_scratch);
  return log_value_;
}

void SlaterDet::mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const RefVector<ParticleSet::ParticleGradient>& G_list,
                              const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                              bool fromscratch) const
{
  constexpr LogValue czero(0);

  for (int iw = 0; iw < wfc_list.size(); iw++)
    wfc_list.getCastedElement<SlaterDet>(iw).log_value_ = czero;

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto Det_list(extract_DetRef_list(wfc_list, i));
    Dets[i]->mw_evaluateGL(Det_list, p_list, G_list, L_list, fromscratch);
    for (int iw = 0; iw < wfc_list.size(); iw++)
      wfc_list.getCastedElement<SlaterDet>(iw).log_value_ += Det_list[iw].get_log_value();
  }
}

void SlaterDet::recompute(const ParticleSet& P)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->recompute(P);
}

void SlaterDet::mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             const std::vector<bool>& recompute) const
{
  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto Det_list(extract_DetRef_list(wfc_list, i));
    Dets[i]->mw_recompute(Det_list, p_list, recompute);
  }
}

void SlaterDet::evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi)
{
  grad_grad_psi.resize(P.getTotalNum());
  HessVector tmp;
  tmp.resize(P.getTotalNum());
  for (int i = 0; i < Dets.size(); ++i)
  {
    tmp = 0;
    Dets[i]->evaluateHessian(P, tmp);
    //  app_log()<<"squee ----- "<<i<< std::endl;
    //  app_log()<<"grad_grad_psi = "<<grad_grad_psi<< std::endl;
    //  app_log()<<"tmp = "<<tmp<< std::endl<< std::endl;
    grad_grad_psi += tmp;
  }
}

void SlaterDet::createResource(ResourceCollection& collection) const
{
  for (const auto& sposet : sposets_)
    sposet->createResource(collection);
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->createResource(collection);
}

void SlaterDet::acquireResource(ResourceCollection& collection,
                                const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& sd_leader = wfc_list.getCastedLeader<SlaterDet>();

  for (int i = 0; i < sposets_.size(); i++)
  {
    auto& phi_leader = *sd_leader.sposets_[i];
    RefVectorWithLeader<SPOSet> phi_list(phi_leader);
    for (WaveFunctionComponent& wfc : wfc_list)
    {
      auto& sd = static_cast<SlaterDet&>(wfc);
      phi_list.push_back(*sd.sposets_[i]);
    }
    phi_leader.acquireResource(collection, phi_list);
  }

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto det_list(extract_DetRef_list(wfc_list, i));
    Dets[i]->acquireResource(collection, det_list);
  }
}

void SlaterDet::releaseResource(ResourceCollection& collection,
                                const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& sd_leader = wfc_list.getCastedLeader<SlaterDet>();

  for (int i = 0; i < sposets_.size(); i++)
  {
    auto& phi_leader = *sd_leader.sposets_[i];
    RefVectorWithLeader<SPOSet> phi_list(phi_leader);
    for (WaveFunctionComponent& wfc : wfc_list)
    {
      auto& sd = static_cast<SlaterDet&>(wfc);
      phi_list.push_back(*sd.sposets_[i]);
    }

    phi_leader.releaseResource(collection, phi_list);
  }

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto det_list(extract_DetRef_list(wfc_list, i));
    Dets[i]->releaseResource(collection, det_list);
  }
}

void SlaterDet::registerData(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
}

SlaterDet::LogValue SlaterDet::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  log_value_ = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    log_value_ += Dets[i]->updateBuffer(P, buf, fromscratch);
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  return log_value_;
}

void SlaterDet::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
}

std::unique_ptr<WaveFunctionComponent> SlaterDet::makeClone(ParticleSet& tqp) const
{
  std::vector<std::unique_ptr<SPOSet>> sposet_clones;
  for (const auto& phi : sposets_)
    sposet_clones.emplace_back(phi->makeClone());

  std::vector<std::unique_ptr<Determinant_t>> dets;
  for (const auto& det : Dets)
  {
    auto it = std::find_if(sposets_.begin(), sposets_.end(),
                           [&](const std::unique_ptr<SPOSet>& sposet) { return sposet.get() == &det->getPhi(); });
    if (it == sposets_.end())
      throw std::runtime_error("Bug! The sposet of a determinant doesn't reference sposets owned by SlaterDet.");
    else
      dets.emplace_back(det->makeCopy(**(sposet_clones.begin() + std::distance(sposets_.begin(), it))));
  }
  auto myclone = std::make_unique<SlaterDet>(tqp, std::move(sposet_clones), std::move(dets));
  assert(myclone->isOptimizable() == isOptimizable());
  return myclone;
}

void SlaterDet::registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const
{
  for (int i = 0; i < Dets.size(); ++i)
  {
    Dets[i]->registerTWFFastDerivWrapper(P, twf);
  }
}

} // namespace qmcplusplus
