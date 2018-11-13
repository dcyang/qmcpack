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
    
    
#ifndef QMCPLUSPLUS_GAUSSIAN_JASTROW_BUILDER_H
#define QMCPLUSPLUS_GAUSSIAN_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{

/** JastrowBuilder using Gaussian functor
 *
 * To replace GaussianConstraints
 */
struct GaussianJastrowBuilder: public WaveFunctionComponentBuilder
{
  ParticleSet *sourcePtcl;
  GaussianJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet& source);
  bool put(xmlNodePtr cur);
  // PtclPoolType& ptclPool;
};
}

#endif
