//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: ChangMo Yang, dcyang@unist.ac.kr, UNIST, Korea
//
// File created by: ChangMo Yang, dcyang@unist.ac.kr, UNIST, Korea
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/GaussianJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/GaussianFunctors.h"
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"

namespace qmcplusplus
{

  GaussianJastrowBuilder::GaussianJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet& source):
    WaveFunctionComponentBuilder(target,psi), sourcePtcl(&source)
  {
    ClassName="GaussianJastrowBuilder";
  }
  /*
  GaussianJastrowBuilder::GaussianJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
      PtclPoolType& psets):
    WaveFunctionComponentBuilder(target,psi),ptclPool(psets)
  {
    ClassName="GaussianJastrowBuilder";
  }
  */

  bool GaussianJastrowBuilder::put(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"put()");
    std::string sourceOpt=targetPtcl.getName();
    std::string jname="GaussianJastrow";
    std::string spin="yes";
    std::string id_a="jab_a";
    OhmmsAttributeSet pattrib;
    pattrib.add(jname,"name");
    pattrib.add(spin,"spin");
    pattrib.add(sourceOpt,"source");
    pattrib.put(cur);
    //bool spindep=(spin=="yes");
    SpeciesSet &sSet = sourcePtcl->getSpeciesSet();
    SpeciesSet &tSet = targetPtcl.getSpeciesSet();
    // int chargeInd=sSet.addAttribute("charge");
    typedef GaussianFunctor<RealType> RadFuncType;
    if(sourceOpt == sourcePtcl->getName())
    {
      //one-body
#if defined(ENABLE_SOA)
      typedef J1OrbitalSoA<RadFuncType> J1Type;
#else
      typedef OneBodyJastrowOrbital<RadFuncType> J1Type;
#endif
      typedef DiffOneBodyJastrowOrbital<RadFuncType> dJ1Type;
      // int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
      J1Type *J1 = new J1Type(*sourcePtcl, targetPtcl);
      dJ1Type *dJ1 = new dJ1Type(*sourcePtcl, targetPtcl);
      cur= cur->xmlChildrenNode;
      while(cur!=NULL)
      {
        std::string cname((const char*)cur->name);
        if (cname == "correlation")
        {
          OhmmsAttributeSet rAttrib;
          RealType cusp=-1e10;
          std::string spA(sSet.speciesName[0]);
          std::string spB(tSet.speciesName[0]);
          rAttrib.add(spA,"speciesA");
          rAttrib.add(spB,"speciesB");
          rAttrib.add(cusp,"cusp");
          rAttrib.put(cur);
          int ia = sSet.findSpecies(spA);
          int ib = tSet.findSpecies(spB);
          if(ia==sSet.size() || ib == tSet.size())
          {
            PRE.error("Failed. Species are incorrect.",true);
          }
          /*
          // TODO: implement as necessary
          if(cusp<-1e6)
          {
            RealType qq=species(chargeInd,ia)*species(chargeInd,ib);
            cusp = (ia==ib)? -0.25*qq:-0.5*qq;
          }
          */
          // std::ostringstream o; o<<"j2"<<ia<<ib;     // this line seems unnecessary
          RadFuncType *functor = new RadFuncType();
          functor->cutoff_radius = targetPtcl.LRBox.LR_rc;
          functor->put(cur);
          J1->addFunc(ia,functor,ib);
          dJ1->addFunc(ia,functor,ib);
          //{
          //  std::ostringstream o;
          //  o<< "gausslan"<<ia<<"-"<<ib;
          //  std::ofstream fstream(o.str().c_str());
          //  int n=100;
          //  RealType d=10/100.,r=0.001;
          //  RealType u,du,d2du;
          //  for (int i=0; i<n; ++i)
          //  {
          //    u=functor->evaluate(r,du,d2du);
          //    fstream << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du
          //      << std::setw(22) << d2du << std::endl;
          //    r+=d;
          //  }
          //}
        }
        cur=cur->next;
      }
      J1->dPsi=dJ1;
      targetPsi.addOrbital(J1,"J1_gaussian");
      J1->setOptimizable(true);
    }
    return true;
  }
}
