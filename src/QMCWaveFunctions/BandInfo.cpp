//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BandInfo.cpp
 */
#include "BandInfo.h"
#include "QMCWaveFunctions/SPOSetInfo.h"
namespace qmcplusplus
{
BandInfoGroup::BandInfoGroup() : NumSPOs(0), FirstBand(0) {}

void BandInfoGroup::selectBands(const std::vector<BandInfo>& bigspace, int first_orb, int num_spos)
{
  app_log() << "BandInfoGroup::selectBands bigspace has " << bigspace.size() << " distinct orbitals " << std::endl;
  myBands.clear();

  int iorb    = 0;
  int N       = bigspace.size();
  int n_lower = 0;
  do
  {
    if (iorb >= first_orb)
      break;
    n_lower += (bigspace[iorb].MakeTwoCopies) ? 2 : 1;
    ++iorb;
  } while (iorb < N);

  if (iorb >= N)
  {
    APP_ABORT("BandInfoGroup::selectBands failed due to iorb>=N");
  }

  FirstBand  = iorb;
  NumSPOs    = 0;
  int ns_max = num_spos - 1;
  while (iorb < N && NumSPOs < num_spos)
  {
    //if(myBands.size()>=num_orbs) break;
    myBands.push_back(bigspace[iorb]);
    NumSPOs += (NumSPOs < ns_max && bigspace[iorb].MakeTwoCopies) ? 2 : 1;
    ++iorb;
  }

  app_log() << "BandInfoGroup::selectBands using distinct orbitals [" << first_orb << "," << iorb << ")" << std::endl;
  app_log() << "  Number of distinct bands " << myBands.size() << std::endl;
  app_log() << "  First Band index " << FirstBand << std::endl;
  app_log() << "  Size of SPOs " << NumSPOs << std::endl;
}

} // namespace qmcplusplus
