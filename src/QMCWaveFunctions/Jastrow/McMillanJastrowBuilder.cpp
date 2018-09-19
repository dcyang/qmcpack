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
    
    
#include "QMCWaveFunctions/Jastrow/McMillanJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"

namespace qmcplusplus
{

  McMillanJastrowBuilder::McMillanJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
      PtclPoolType& psets):
    WaveFunctionComponentBuilder(target,psi),ptclPool(psets)
  {
    ClassName="McMillanJastrowBuilder";
  }

  bool McMillanJastrowBuilder::put(xmlNodePtr cur)
  {
    ReportEngine PRE(ClassName,"put()");
    std::string sourceOpt=targetPtcl.getName();
    std::string jname="McMillanJastrow";
    std::string spin="yes";
    std::string id_a="jaa_a";
    std::string id_b="jaa_b";
    RealType pade_a=1.0;
    RealType pade_b=1.0;
    OhmmsAttributeSet pattrib;
    pattrib.add(jname,"name");
    pattrib.add(spin,"spin");
    pattrib.add(sourceOpt,"source");
    pattrib.put(cur);
    //bool spindep=(spin=="yes");
    SpeciesSet& species(targetPtcl.getSpeciesSet());
    int chargeInd=species.addAttribute("charge");
    typedef McMillanFunctor<RealType> RadFuncType;
    if(sourceOpt == targetPtcl.getName())
    {
      //two-body
#if defined(ENABLE_SOA)
      typedef J2OrbitalSoA<RadFuncType> J2Type;
#else
      typedef TwoBodyJastrowOrbital<RadFuncType> J2Type;
#endif
      typedef DiffTwoBodyJastrowOrbital<RadFuncType> dJ2Type;
      int taskid=(targetPsi.is_manager())?targetPsi.getGroupID():-1;
      J2Type *J2 = new J2Type(targetPtcl,taskid);
      dJ2Type *dJ2 = new dJ2Type(targetPtcl);
      cur= cur->xmlChildrenNode;
      while(cur!=NULL)
      {
        std::string cname((const char*)cur->name);
        if (cname == "correlation")
        {
          OhmmsAttributeSet rAttrib;
          RealType cusp=-1e10;
          std::string spA(species.speciesName[0]);
          std::string spB(species.speciesName[0]);
          rAttrib.add(spA,"speciesA");
          rAttrib.add(spB,"speciesB");
          rAttrib.add(cusp,"cusp");
          rAttrib.put(cur);
          int ia = species.findSpecies(spA);
          int ib = species.findSpecies(spB);
          if(ia==species.size() || ib == species.size())
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
          std::ostringstream o;
          o<<"j2"<<ia<<ib;
          RadFuncType *functor = new RadFuncType();
          functor->cutoff_radius = targetPtcl.LRBox.LR_rc;
          functor->put(cur);
          J2->addFunc(ia,ib,functor);
          dJ2->addFunc(ia,ib,functor);
          //{
          //  std::ostringstream o;
          //  o<< "mcmillan"<<ia<<"-"<<ib;
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
      J2->dPsi=dJ2;
      targetPsi.addOrbital(J2,"J2_mcmillan");
      J2->setOptimizable(true);
    }
    return true;
  }
}
