//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCCostFunction.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
//#define QMCCOSTFUNCTION_DEBUG

namespace qmcplusplus
{
QMCCostFunction::QMCCostFunction(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm)
    : QMCCostFunctionBase(w, psi, h, comm),
      fill_timer_(createGlobalTimer("QMCCostFunction::fillOverlapHamiltonianMatrices", timer_level_medium))
{
  CSWeight = 1.0;
  app_log() << " Using QMCCostFunction::QMCCostFunction" << std::endl;
}


/** Clean up the vector */
QMCCostFunction::~QMCCostFunction()
{
  delete_iter(RecordsOnNode.begin(), RecordsOnNode.end());
  delete_iter(DerivRecords.begin(), DerivRecords.end());
  delete_iter(HDerivRecords.begin(), HDerivRecords.end());
}

void QMCCostFunction::GradCost(std::vector<Return_rt>& PGradient,
                               const std::vector<Return_rt>& PM,
                               Return_rt FiniteDiff)
{
  for (int j = 0; j < NumOptimizables; j++)
    OptVariables[j] = PM[j];
  if (FiniteDiff > 0)
  {
    QMCTraits::RealType dh = 1.0 / (2.0 * FiniteDiff);
    for (int i = 0; i < NumOptimizables; i++)
    {
      // + FiniteDiff
      OptVariables[i] = PM[i] + FiniteDiff;
      resetPsi();
      correlatedSampling(false);
      auto CostPlus = computedCost();
      // - FiniteDiff
      OptVariables[i] = PM[i] - FiniteDiff;
      resetPsi();
      correlatedSampling(false);
      auto CostMinus = computedCost();
      // calculate gradient
      PGradient[i]    = (CostPlus - CostMinus) * dh;
      OptVariables[i] = PM[i]; // revert parameter change
    }
  }
  else
  {
    resetPsi();
    //evaluate new local energies and derivatives
    EffectiveWeight effective_weight = correlatedSampling(true);
    //Estimators::accumulate has been called by correlatedSampling
    curAvg_w = SumValue[SUM_E_WGT] / SumValue[SUM_WGT];
    //    Return_t curAvg2_w = curAvg_w*curAvg_w;
    curVar_w = SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT] - curAvg_w * curAvg_w;
    std::vector<Return_rt> EDtotals(NumOptimizables, 0.0);
    std::vector<Return_rt> EDtotals_w(NumOptimizables, 0.0);
    std::vector<Return_rt> E2Dtotals_w(NumOptimizables, 0.0);
    std::vector<Return_rt> URV(NumOptimizables, 0.0);
    std::vector<Return_rt> HD_avg(NumOptimizables, 0.0);
    Return_rt wgtinv   = 1.0 / SumValue[SUM_WGT];
    Return_rt delE_bar = 0;
    for (int ip = 0; ip < NumThreads; ip++)
    {
      int nw = wClones[ip]->numSamples();
      for (int iw = 0; iw < nw; iw++)
      {
        const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        delE_bar += weight * std::pow(std::abs(eloc_new - EtargetEff), PowerE);
        const Return_rt* HDsaved = (*HDerivRecords[ip])[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
          HD_avg[pm] += HDsaved[pm];
      }
    }
    myComm->allreduce(HD_avg);
    myComm->allreduce(delE_bar);
    for (int pm = 0; pm < NumOptimizables; pm++)
      HD_avg[pm] *= 1.0 / static_cast<Return_rt>(NumSamples);
    for (int ip = 0; ip < NumThreads; ip++)
    {
      int nw = wClones[ip]->numSamples();
      for (int iw = 0; iw < nw; iw++)
      {
        const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        Return_rt delta_l               = (eloc_new - curAvg_w);
        bool ltz(true);
        if (eloc_new - EtargetEff < 0)
          ltz = false;
        Return_rt delE           = std::pow(std::abs(eloc_new - EtargetEff), PowerE);
        Return_rt ddelE          = PowerE * std::pow(std::abs(eloc_new - EtargetEff), PowerE - 1);
        const Return_t* Dsaved   = (*DerivRecords[ip])[iw];
        const Return_rt* HDsaved = (*HDerivRecords[ip])[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
        {
          EDtotals_w[pm] += weight * (HDsaved[pm] + 2.0 * std::real(Dsaved[pm]) * delta_l);
          URV[pm] += 2.0 * (eloc_new * HDsaved[pm] - curAvg * HD_avg[pm]);
          if (ltz)
            EDtotals[pm] += weight * (2.0 * std::real(Dsaved[pm]) * (delE - delE_bar) + ddelE * HDsaved[pm]);
          else
            EDtotals[pm] += weight * (2.0 * std::real(Dsaved[pm]) * (delE - delE_bar) - ddelE * HDsaved[pm]);
        }
      }
    }
    myComm->allreduce(EDtotals);
    myComm->allreduce(EDtotals_w);
    myComm->allreduce(URV);
    Return_rt smpinv = 1.0 / static_cast<Return_rt>(NumSamples);
    for (int ip = 0; ip < NumThreads; ip++)
    {
      int nw = wClones[ip]->numSamples();
      for (int iw = 0; iw < nw; iw++)
      {
        const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_rt weight                = saved[REWEIGHT] * wgtinv;
        Return_rt eloc_new              = saved[ENERGY_NEW];
        Return_rt delta_l               = (eloc_new - curAvg_w);
        Return_rt sigma_l               = delta_l * delta_l;
        const Return_t* Dsaved          = (*DerivRecords[ip])[iw];
        const Return_rt* HDsaved        = (*HDerivRecords[ip])[iw];
        for (int pm = 0; pm < NumOptimizables; pm++)
        {
          E2Dtotals_w[pm] +=
              weight * 2.0 * (std::real(Dsaved[pm]) * (sigma_l - curVar_w) + delta_l * (HDsaved[pm] - EDtotals_w[pm]));
        }
      }
    }
    myComm->allreduce(E2Dtotals_w);
    for (int pm = 0; pm < NumOptimizables; pm++)
      URV[pm] *= smpinv;
    for (int j = 0; j < NumOptimizables; j++)
    {
      PGradient[j] = 0.0;
      if (std::abs(w_var) > 1.0e-10)
        PGradient[j] += w_var * E2Dtotals_w[j];
      if (std::abs(w_en) > 1.0e-10)
        PGradient[j] += w_en * EDtotals_w[j];
      if (std::abs(w_w) > 1.0e-10)
        PGradient[j] += w_w * URV[j];
      if (std::abs(w_abs) > 1.0e-10)
        PGradient[j] += w_abs * EDtotals[j];
    }

    IsValid = isEffectiveWeightValid(effective_weight);
  }
}


void QMCCostFunction::getConfigurations(const std::string& aroot)
{
  //makeClones(W,Psi,H);
  if (H_KE_Node.empty())
  {
    app_log() << "  QMCCostFunction is created with " << NumThreads << " threads." << std::endl;
    //make H_KE_Node
    H_KE_Node.resize(NumThreads);
    RecordsOnNode.resize(NumThreads, 0);
    DerivRecords.resize(NumThreads, 0);
    HDerivRecords.resize(NumThreads, 0);
  }

  //#pragma omp parallel for
  for (int ip = 0; ip < NumThreads; ++ip)
    if (!H_KE_Node[ip])
    {
      auto components = hClones[ip]->getTWFDependentComponents();
      if (ip == 0)
      {
        app_log() << " Found " << components.size() << " wavefunction dependent components in the Hamiltonian";
        if (components.size())
          for (const OperatorBase& component : components)
            app_log() << " '" << component.getName() << "'";
        app_log() << "." << std::endl;
      }
      H_KE_Node[ip] = std::make_unique<HamiltonianRef>(components);
    }

  //load samples from SampleStack
  app_log() << "   Number of samples loaded to each thread : ";
  wPerRank[0] = 0;
  for (int ip = 0; ip < NumThreads; ++ip)
  {
    wPerRank[ip + 1] = wPerRank[ip] + wClones[ip]->numSamples();
    app_log() << wClones[ip]->numSamples() << " ";
  }
  app_log() << std::endl;
  app_log().flush();

  if (dLogPsi.size() != wPerRank[NumThreads])
  {
    delete_iter(dLogPsi.begin(), dLogPsi.end());
    delete_iter(d2LogPsi.begin(), d2LogPsi.end());
    int nptcl = W.getTotalNum();
    int nwtot = wPerRank[NumThreads];
    dLogPsi.resize(nwtot);
    d2LogPsi.resize(nwtot);
    for (int i = 0; i < nwtot; ++i)
      dLogPsi[i] = new ParticleGradient(nptcl);
    for (int i = 0; i < nwtot; ++i)
      d2LogPsi[i] = new ParticleLaplacian(nptcl);
  }
}

/** evaluate everything before optimization */
void QMCCostFunction::checkConfigurations(EngineHandle& handle)
{
  RealType et_tot = 0.0;
  RealType e2_tot = 0.0;
#pragma omp parallel reduction(+ : et_tot, e2_tot)
  {
    int ip = omp_get_thread_num();
    MCWalkerConfiguration& wRef(*wClones[ip]);
    if (RecordsOnNode[ip] == 0)
    {
      RecordsOnNode[ip] = new Matrix<Return_rt>;
      RecordsOnNode[ip]->resize(wRef.numSamples(), SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip] = new Matrix<Return_t>;
        DerivRecords[ip]->resize(wRef.numSamples(), NumOptimizables);
        HDerivRecords[ip] = new Matrix<Return_rt>;
        HDerivRecords[ip]->resize(wRef.numSamples(), NumOptimizables);
      }
    }
    else if (RecordsOnNode[ip]->size1() != wRef.numSamples())
    {
      RecordsOnNode[ip]->resize(wRef.numSamples(), SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]->resize(wRef.numSamples(), NumOptimizables);
        HDerivRecords[ip]->resize(wRef.numSamples(), NumOptimizables);
      }
    }
    // Populate local to global index mapping into psiClone internal component 'myVars',
    // because psiClones persist between different sections and need update.
    psiClones[ip]->checkOutVariables(OptVariablesForPsi);
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
    //int nat = wRef.getTotalNum();
    Return_rt e0 = 0.0;
    //       Return_t ef=0.0;
    Return_rt e2 = 0.0;
    for (int iw = 0, iwg = wPerRank[ip]; iw < wRef.numSamples(); ++iw, ++iwg)
    {
      wRef.loadSample(wRef, iw);
      wRef.update();
      Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      psiClones[ip]->evaluateDeltaLogSetup(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg],
                                           *d2LogPsi[iwg]);
      saved[REWEIGHT] = 1.0;
      Return_rt etmp;
      if (needGrads)
      {
        //allocate vector
        Vector<Return_rt> rDsaved(NumOptimizables, 0.0);
        Vector<Return_rt> rHDsaved(NumOptimizables, 0.0);

        Vector<Return_t> Dsaved(NumOptimizables, 0.0);
        Vector<Return_t> HDsaved(NumOptimizables, 0.0);

        etmp = hClones[ip]->evaluateValueAndDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);


        //FIXME the ifdef should be removed after the optimizer is made compatible with complex coefficients
        for (int i = 0; i < NumOptimizables; i++)
        {
          rHDsaved[i] = std::real(HDsaved[i]);
        }
        std::copy(Dsaved.begin(), Dsaved.end(), (*DerivRecords[ip])[iw]);
        std::copy(rHDsaved.begin(), rHDsaved.end(), (*HDerivRecords[ip])[iw]);
      }
      else
        etmp = hClones[ip]->evaluate(wRef);

      e0 += saved[ENERGY_TOT] = saved[ENERGY_NEW] = etmp;
      e2 += etmp * etmp;
      saved[ENERGY_FIXED]                 = saved[ENERGY_TOT];
      const auto twf_dependent_components = hClones[ip]->getTWFDependentComponents();
      for (const OperatorBase& component : twf_dependent_components)
        saved[ENERGY_FIXED] -= component.getValue();
    }
    //add them all using reduction
    et_tot += e0;
    e2_tot += e2;
    // #pragma omp atomic
    //       eft_tot+=ef;
  }
  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerRank[NumThreads]) << std::endl;
  //Need to sum over the processors
  std::vector<Return_rt> etemp(3);
  etemp[0] = et_tot;
  etemp[1] = static_cast<Return_rt>(wPerRank[NumThreads]);
  etemp[2] = e2_tot;
  myComm->allreduce(etemp);
  Etarget    = static_cast<Return_rt>(etemp[0] / etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << std::endl;
  app_log() << "  VMC Evar = " << etemp[2] / etemp[1] - Etarget * Etarget << std::endl;
  app_log() << "  Total weights = " << etemp[1] << std::endl;
  app_log().flush();
  setTargetEnergy(Etarget);
  ReportCounter = 0;
  IsValid       = true;
  //collect SumValue for computedCost
  SumValue[SUM_WGT]       = etemp[1];
  SumValue[SUM_WGTSQ]     = etemp[1];
  SumValue[SUM_E_WGT]     = etemp[0];
  SumValue[SUM_ESQ_WGT]   = etemp[2];
  SumValue[SUM_E_BARE]    = etemp[0];
  SumValue[SUM_ESQ_BARE]  = etemp[2];
  SumValue[SUM_ABSE_BARE] = 0.0;
}

#ifdef HAVE_LMY_ENGINE
/** evaluate everything before optimization
 *In future, both the LM and descent engines should be children of some parent engine base class.
 * */
void QMCCostFunction::engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                                 DescentEngine& descentEngineObj,
                                                 const std::string& MinMethod)
{
  if (MinMethod == "descent")
  {
    //Reset vectors and scalars from any previous iteration
    descentEngineObj.prepareStorage(omp_get_max_threads(), NumOptimizables);
  }
  RealType et_tot = 0.0;
  RealType e2_tot = 0.0;
#pragma omp parallel reduction(+ : et_tot, e2_tot)
  {
    int ip = omp_get_thread_num();
    MCWalkerConfiguration& wRef(*wClones[ip]);
    if (RecordsOnNode[ip] == 0)
    {
      RecordsOnNode[ip] = new Matrix<Return_rt>;
      RecordsOnNode[ip]->resize(wRef.numSamples(), SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip] = new Matrix<Return_t>;
        //DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        HDerivRecords[ip] = new Matrix<Return_rt>;
        //HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    else if (RecordsOnNode[ip]->size1() != wRef.numSamples())
    {
      RecordsOnNode[ip]->resize(wRef.numSamples(), SUM_INDEX_SIZE);
      if (needGrads)
      {
        //DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        //HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    // Populate local to global index mapping into psiClone internal component 'myVars',
    // because psiClones persist between different sections and need update.
    psiClones[ip]->checkOutVariables(OptVariablesForPsi);
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
    //int nat = wRef.getTotalNum();
    Return_rt e0 = 0.0;
    //       Return_t ef=0.0;
    Return_rt e2 = 0.0;


    for (int iw = 0, iwg = wPerRank[ip]; iw < wRef.numSamples(); ++iw, ++iwg)
    {
      wRef.loadSample(wRef, iw);
      wRef.update();
      Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      psiClones[ip]->evaluateDeltaLogSetup(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg],
                                           *d2LogPsi[iwg]);
      saved[REWEIGHT] = 1.0;
      Return_rt etmp;
      if (needGrads)
      {
        //allocate vector
        Vector<Return_t> Dsaved(NumOptimizables, 0.0);
        Vector<Return_t> HDsaved(NumOptimizables, 0.0);

        etmp = hClones[ip]->evaluateValueAndDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);

        // add non-differentiated derivative vector
        std::vector<Return_t> der_rat_samp(NumOptimizables + 1, 0.0);
        std::vector<Return_t> le_der_samp(NumOptimizables + 1, 0.0);

        // dervative vectors
        der_rat_samp.at(0) = 1.0;
        for (int i = 0; i < Dsaved.size(); i++)
          der_rat_samp[i + 1] = Dsaved[i];

        // energy dervivatives
        le_der_samp.at(0) = etmp;
        for (int i = 0; i < HDsaved.size(); i++)
          le_der_samp[i + 1] = HDsaved[i] + etmp * Dsaved[i];

#ifdef HAVE_LMY_ENGINE
        if (MinMethod == "adaptive")
        {
          // pass into engine
          EngineObj->take_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, saved[REWEIGHT]);
        }
        else if (MinMethod == "descent")
        {
          //Could remove this copying over if LM engine becomes compatible with complex numbers
          //so that der_rat_samp and le_der_samp are vectors of std::complex<double> when QMC_COMPLEX=1
          std::vector<FullPrecValueType> der_rat_samp_comp(der_rat_samp.begin(), der_rat_samp.end());
          std::vector<FullPrecValueType> le_der_samp_comp(le_der_samp.begin(), le_der_samp.end());

          descentEngineObj.takeSample(ip, der_rat_samp_comp, le_der_samp_comp, le_der_samp_comp, 1.0, saved[REWEIGHT]);
        }
#endif
      }
      else
        etmp = hClones[ip]->evaluate(wRef);

      e0 += saved[ENERGY_TOT] = etmp;
      e2 += etmp * etmp;

      saved[ENERGY_FIXED]                 = saved[ENERGY_TOT];
      const auto twf_dependent_components = hClones[ip]->getTWFDependentComponents();
      for (const OperatorBase& component : twf_dependent_components)
        saved[ENERGY_FIXED] -= component.getValue();
    }

    //add them all using reduction
    et_tot += e0;
    e2_tot += e2;
    // #pragma omp atomic
    //       eft_tot+=ef;
  }

  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerRank[NumThreads]) << endl;
  //Need to sum over the processors
  std::vector<Return_rt> etemp(3);
  etemp[0] = et_tot;
  etemp[1] = static_cast<Return_rt>(wPerRank[NumThreads]);
  etemp[2] = e2_tot;
  myComm->allreduce(etemp);
  Etarget    = static_cast<Return_rt>(etemp[0] / etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << std::endl;
  app_log() << "  VMC Evar = " << etemp[2] / etemp[1] - Etarget * Etarget << std::endl;
  app_log() << "  Total weights = " << etemp[1] << std::endl;


#ifdef HAVE_LMY_ENGINE
  // engine finish taking samples
  if (MinMethod == "adaptive")
    EngineObj->sample_finish();
  else if (MinMethod == "descent")
    descentEngineObj.sample_finish();
#endif

  app_log().flush();

  setTargetEnergy(Etarget);
  ReportCounter = 0;
}
#endif


void QMCCostFunction::resetPsi(bool final_reset)
{
  for (int i = 0; i < OptVariables.size(); ++i)
    OptVariablesForPsi[i] = OptVariables[i];
  //cout << "######### QMCCostFunction::resetPsi " << std::endl;
  //OptVariablesForPsi.print(std::cout);
  //cout << "-------------------------------------- " << std::endl;

  resetOptimizableObjects(Psi, OptVariablesForPsi);
  for (int i = 0; i < psiClones.size(); ++i)
    resetOptimizableObjects(*psiClones[i], OptVariablesForPsi);
}

QMCCostFunction::EffectiveWeight QMCCostFunction::correlatedSampling(bool needGrad)
{
  for (int ip = 0; ip < NumThreads; ++ip)
  {
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
  }

  Return_rt wgt_tot       = 0.0;
  Return_rt wgt_tot2      = 0.0;
  Return_rt inv_n_samples = 1.0 / NumSamples;
#pragma omp parallel reduction(+ : wgt_tot, wgt_tot2)
  {
    const int ip = omp_get_thread_num();
    //if we have more than KE depending on TWF, TWF must be fully recomputed.
    const bool compute_all_from_scratch = hClones[ip]->getTWFDependentComponents().size() > 1;

    MCWalkerConfiguration& wRef(*wClones[ip]);
    Return_rt wgt_node = 0.0, wgt_node2 = 0.0;
    for (int iw = 0, iwg = wPerRank[ip]; iw < wRef.numSamples(); ++iw, ++iwg)
    {
      wRef.loadSample(wRef, iw);
      wRef.update(true);
      Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_rt logpsi;
      logpsi = psiClones[ip]->evaluateDeltaLog(wRef, compute_all_from_scratch);
      wRef.G += *dLogPsi[iwg];
      wRef.L += *d2LogPsi[iwg];
      Return_rt weight = saved[REWEIGHT] = vmc_or_dmc * (logpsi - saved[LOGPSI_FREE]);
      if (needGrad)
      {
        Vector<Return_t> Dsaved(NumOptimizables, 0);
        Vector<Return_t> HDsaved(NumOptimizables, 0);

        Vector<Return_rt> rDsaved(NumOptimizables, 0);
        Vector<Return_rt> rHDsaved(NumOptimizables, 0);

        saved[ENERGY_NEW] =
            H_KE_Node[ip]->evaluateValueAndDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved) + saved[ENERGY_FIXED];
        ;

        for (int i = 0; i < NumOptimizables; i++)
        {
          rDsaved[i]  = std::real(Dsaved[i]);
          rHDsaved[i] = std::real(HDsaved[i]);
          (*DerivRecords[ip])(iw, i)  = rDsaved[i];
          (*HDerivRecords[ip])(iw, i) = rHDsaved[i];
        }
      }
      else
        saved[ENERGY_NEW] = H_KE_Node[ip]->evaluate(wRef) + saved[ENERGY_FIXED];
      wgt_node += inv_n_samples * weight;
      wgt_node2 += inv_n_samples * weight * weight;
    }
    wgt_tot += wgt_node;
    wgt_tot2 += wgt_node2;
  }
  //this is MPI barrier
  OHMMS::Controller->barrier();
  //collect the total weight for normalization and apply maximum weight
  myComm->allreduce(wgt_tot);
  myComm->allreduce(wgt_tot2);
  //    app_log()<<"Before Purge"<<wgt_tot<<" "<<wgt_tot2<< std::endl;
  Return_rt wgtnorm = (wgt_tot == 0) ? 0 : wgt_tot;
  wgt_tot           = 0.0;
  for (int ip = 0; ip < NumThreads; ip++)
  {
    int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++)
    {
      Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT] =
          std::min(std::exp(saved[REWEIGHT] - wgtnorm), std::numeric_limits<Return_rt>::max() * (RealType)0.1);
      wgt_tot += inv_n_samples * saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"During Purge"<<wgt_tot<<" "<< std::endl;
  wgtnorm = (wgt_tot == 0) ? 1 : 1.0 / wgt_tot;
  wgt_tot = 0.0;
  for (int ip = 0; ip < NumThreads; ip++)
  {
    int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++)
    {
      Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT]           = std::min(saved[REWEIGHT] * wgtnorm, MaxWeight);
      wgt_tot += inv_n_samples * saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"After Purge"<<wgt_tot<<" "<< std::endl;
  for (int i = 0; i < SumValue.size(); i++)
    SumValue[i] = 0.0;
  CSWeight = wgt_tot = (wgt_tot == 0) ? 1 : 1.0 / wgt_tot;
  for (int ip = 0; ip < NumThreads; ip++)
  {
    int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++)
    {
      const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      //      Return_t weight=saved[REWEIGHT]*wgt_tot;
      Return_rt eloc_new = saved[ENERGY_NEW];
      Return_rt delE     = std::pow(std::abs(eloc_new - EtargetEff), PowerE);
      SumValue[SUM_E_BARE] += eloc_new;
      SumValue[SUM_ESQ_BARE] += eloc_new * eloc_new;
      SumValue[SUM_ABSE_BARE] += delE;
      SumValue[SUM_E_WGT] += eloc_new * saved[REWEIGHT];
      SumValue[SUM_ESQ_WGT] += eloc_new * eloc_new * saved[REWEIGHT];
      SumValue[SUM_ABSE_WGT] += delE * saved[REWEIGHT];
      SumValue[SUM_WGT] += saved[REWEIGHT];
      SumValue[SUM_WGTSQ] += saved[REWEIGHT] * saved[REWEIGHT];
    }
  }
  //collect everything
  myComm->allreduce(SumValue);
  return SumValue[SUM_WGT] * SumValue[SUM_WGT] / (SumValue[SUM_WGTSQ] * NumSamples);
}


QMCCostFunction::Return_rt QMCCostFunction::fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left,
                                                                           Matrix<Return_rt>& Right)
{
  ScopedTimer tmp_timer(fill_timer_);

  RealType b2(w_beta);

  Right = 0.0;
  Left  = 0.0;

  //     resetPsi();
  curAvg_w            = SumValue[SUM_E_WGT] / SumValue[SUM_WGT];
  Return_rt curAvg2_w = SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT];
  RealType V_avg      = curAvg2_w - curAvg_w * curAvg_w;
  std::vector<Return_t> D_avg(getNumParams(), 0.0);
  Return_rt wgtinv = 1.0 / SumValue[SUM_WGT];
  for (int ip = 0; ip < NumThreads; ip++)
  {
    int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++)
    {
      const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_rt weight                = saved[REWEIGHT] * wgtinv;
      const Return_t* Dsaved          = (*DerivRecords[ip])[iw];
      for (int pm = 0; pm < getNumParams(); pm++)
      {
        D_avg[pm] += Dsaved[pm] * weight;
      }
    }
  }

  myComm->allreduce(D_avg);

  for (int ip = 0; ip < NumThreads; ip++)
  {
    int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++)
    {
      const Return_rt* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_rt weight                = saved[REWEIGHT] * wgtinv;
      Return_rt eloc_new              = saved[ENERGY_NEW];
      const Return_t* Dsaved          = (*DerivRecords[ip])[iw];
      const Return_rt* HDsaved        = (*HDerivRecords[ip])[iw];
#pragma omp parallel for
      for (int pm = 0; pm < getNumParams(); pm++)
      {
        Return_t wfe   = (HDsaved[pm] + (Dsaved[pm] - D_avg[pm]) * eloc_new) * weight;
        Return_t wfd   = (Dsaved[pm] - D_avg[pm]) * weight;
        Return_t vterm = HDsaved[pm] * (eloc_new - curAvg_w) +
            (Dsaved[pm] - D_avg[pm]) * eloc_new * (eloc_new - RealType(2.0) * curAvg_w);
        //                Return_t vterm = (HDsaved[pm]+(Dsaved[pm]-D_avg[pm])*eloc_new -curAvg_w)*(eloc_new-curAvg_w);
        //                 Variance
        Left(0, pm + 1) += b2 * std::real(vterm) * weight;
        Left(pm + 1, 0) += b2 * std::real(vterm) * weight;
        //                 Hamiltonian
        Left(0, pm + 1) += (1 - b2) * std::real(wfe);
        Left(pm + 1, 0) += (1 - b2) * std::real(wfd) * eloc_new;
        for (int pm2 = 0; pm2 < getNumParams(); pm2++)
        {
          //                Hamiltonian
          Left(pm + 1, pm2 + 1) +=
              std::real((1 - b2) * std::conj(wfd) * (HDsaved[pm2] + (Dsaved[pm2] - D_avg[pm2]) * eloc_new));
          //                Overlap
          RealType ovlij = std::real(std::conj(wfd) * (Dsaved[pm2] - D_avg[pm2]));
          Right(pm + 1, pm2 + 1) += ovlij;
          //                Variance
          RealType varij = weight *
              std::real((HDsaved[pm] - RealType(2.0) * std::conj(Dsaved[pm] - D_avg[pm]) * eloc_new) *
                        (HDsaved[pm2] - RealType(2.0) * (Dsaved[pm2] - D_avg[pm2]) * eloc_new));
          //                  RealType varij=weight*(HDsaved[pm] +(Dsaved[pm]-D_avg[pm])*eloc_new-curAvg_w)*
          //                                      (HDsaved[pm2] + (Dsaved[pm2]-D_avg[pm2])*eloc_new-curAvg_w);
          Left(pm + 1, pm2 + 1) += b2 * (varij + V_avg * ovlij);
        }
      }
    }
  }
  myComm->allreduce(Right);
  myComm->allreduce(Left);
  Left(0, 0)  = (1 - b2) * curAvg_w + b2 * V_avg;
  Right(0, 0) = 1.0;

  return 1.0;
}
} // namespace qmcplusplus
