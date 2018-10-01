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
    
    
#ifndef QMCPLUSPLUS_MCMILLAN_JASTROW_BUILDER_H
#define QMCPLUSPLUS_MCMILLAN_JASTROW_BUILDER_H
#include "QMCWaveFunctions/Jastrow/McMillanFunctors.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{

/** JastrowBuilder using McMillan functor
 *
 * To replace McMillanConstraints
 */
struct McMillanJastrowBuilder: public WaveFunctionComponentBuilder
{

  McMillanJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, PtclPoolType& psets);
  bool put(xmlNodePtr cur);
  PtclPoolType& ptclPool;
};
}

#endif
