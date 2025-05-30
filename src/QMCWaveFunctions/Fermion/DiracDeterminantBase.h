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


/**@file DiracDeterminantBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_BASE_H
#define QMCPLUSPLUS_DIRACDETERMINANT_BASE_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/TimerManager.h"
#include "CPU/math.hpp"

namespace qmcplusplus
{
/// determinant matrix inverter select
enum class DetMatInvertor
{
  HOST,
  ACCEL,
};

class DiracDeterminantBase : public WaveFunctionComponent
{
public:
  /** constructor
   *@param spos the single-particle orbital set.
   *@param first index of the first particle
   *@param last index of last particle
   */
  DiracDeterminantBase(const std::string& class_name, SPOSet& phi, int first, int last)
      : UpdateTimer(createGlobalTimer(class_name + "::update", timer_level_fine)),
        RatioTimer(createGlobalTimer(class_name + "::ratio", timer_level_fine)),
        InverseTimer(createGlobalTimer(class_name + "::inverse", timer_level_fine)),
        BufferTimer(createGlobalTimer(class_name + "::buffer", timer_level_fine)),
        SPOVTimer(createGlobalTimer(class_name + "::spoval", timer_level_fine)),
        SPOVGLTimer(createGlobalTimer(class_name + "::spovgl", timer_level_fine)),
        phi_(phi),
        FirstIndex(first),
        LastIndex(last),
        NumOrbitals(last - first),
        NumPtcls(last - first)
  {}

  ///default destructor
  ~DiracDeterminantBase() override {}

  // copy constructor and assign operator disabled
  DiracDeterminantBase(const DiracDeterminantBase& s)            = delete;
  DiracDeterminantBase& operator=(const DiracDeterminantBase& s) = delete;

  // get the SPO pointer
  inline SPOSet& getPhi() { return phi_; }

  // get FirstIndex, Last Index
  inline int getFirstIndex() const { return FirstIndex; }
  inline int getLastIndex() const { return LastIndex; }

#ifndef NDEBUG
  virtual ValueMatrix& getPsiMinv() { return dummy_vmt; }
#endif

  bool isFermionic() const final { return true; }
  inline bool isOptimizable() const final { return phi_.isOptimizable(); }

  virtual void registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const override
  {
    throw std::runtime_error("DiracDeterminantBase::registerTWFFastDerivWrapper must be overridden\n");
  }

  virtual void evaluateDerivativesWF(ParticleSet& P,
                                     const opt_variables_type& optvars,
                                     Vector<ValueType>& dlogpsi) override
  {
    // assume no orbital optimization. If implemented, override this function
  }

  // expose CPU interfaces
  using WaveFunctionComponent::evaluateDerivatives;
  using WaveFunctionComponent::evaluateGL;
  using WaveFunctionComponent::evaluateLog;
  using WaveFunctionComponent::mw_evaluateGL;
  using WaveFunctionComponent::mw_evaluateLog;
  using WaveFunctionComponent::recompute;

  using WaveFunctionComponent::copyFromBuffer;
  using WaveFunctionComponent::registerData;
  using WaveFunctionComponent::updateBuffer;

  using WaveFunctionComponent::acquireResource;
  using WaveFunctionComponent::createResource;
  using WaveFunctionComponent::releaseResource;

  using WaveFunctionComponent::acceptMove;
  using WaveFunctionComponent::completeUpdates;
  using WaveFunctionComponent::evalGrad;
  using WaveFunctionComponent::mw_accept_rejectMove;
  using WaveFunctionComponent::mw_calcRatio;
  using WaveFunctionComponent::mw_completeUpdates;
  using WaveFunctionComponent::mw_evalGrad;
  using WaveFunctionComponent::mw_ratioGrad;
  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::ratioGrad;
  using WaveFunctionComponent::restore;

  using WaveFunctionComponent::evalGradSource;
  using WaveFunctionComponent::evaluateHessian;
  using WaveFunctionComponent::evaluateRatios;
  using WaveFunctionComponent::evaluateRatiosAlltoOne;
  using WaveFunctionComponent::evaluateSpinorRatios;
  using WaveFunctionComponent::mw_evaluateRatios;

  inline virtual void mw_evaluateSpinorRatios(
      const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
      const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
      const RefVector<std::pair<ValueVector, ValueVector>>& spinor_multiplier_list,
      std::vector<std::vector<ValueType>>& ratios) const override
  {
    mw_evaluateSpinorRatios_serialized(wfc_list, vp_list, spinor_multiplier_list, ratios);
  }

  // used by DiracDeterminantWithBackflow
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   int offset,
                                   Matrix<RealType>& dlogpsi,
                                   Array<GradType, 3>& dG,
                                   Matrix<RealType>& dL)
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::evaluateDerivatives");
  }

  // final keyword is intended to disable makeClone being further inherited.
  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const final
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::makeClone");
    return std::unique_ptr<DiracDeterminantBase>();
  }

  PsiValue ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad) override
  {
    APP_ABORT("  DiracDeterminantBase::ratioGradWithSpin():  Implementation required\n");
    return 0.0;
  }
  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) override
  {
    APP_ABORT("  DiracDeterminantBase::evalGradWithSpin():  Implementation required\n");
    return GradType();
  }
  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  virtual std::unique_ptr<DiracDeterminantBase> makeCopy(SPOSet& phi) const = 0;

protected:
  /// Timers
  NewTimer &UpdateTimer, &RatioTimer, &InverseTimer, &BufferTimer, &SPOVTimer, &SPOVGLTimer;
  /// a set of single-particle orbitals used to fill in the  values of the matrix
  SPOSet& phi_;
  ///index of the first particle with respect to the particle set
  const int FirstIndex;
  ///index of the last particle with respect to the particle set
  const int LastIndex;
  ///number of single-particle orbitals which belong to this Dirac determinant
  const int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  const int NumPtcls;

#ifndef NDEBUG
  // This is for debugging and testing in debug mode
  // psiMinv is not a base class data member or public in most implementations
  // it is frequently Dual and its consistency not guaranteed.
  ValueMatrix dummy_vmt;
#endif
};

} // namespace qmcplusplus
#endif
