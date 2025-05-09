//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "Lattice/CrystalLattice.h"
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"

namespace qmcplusplus
{
using mRealType = EwaldHandler3D::mRealType;

/** evalaute bare Coulomb using EwaldHandler3D
 */
TEST_CASE("ewald3d", "[lrhandler]")
{
  Lattice lattice;
  lattice.BoxBConds     = true;
  lattice.LR_dim_cutoff = 30.;
  lattice.R.diagonal(5.0);
  lattice.reset();
  CHECK(lattice.Volume == Approx(125));
  lattice.SetLRCutoffs(lattice.Rv);
  //lattice.printCutoffs(app_log());
  CHECK(lattice.LR_rc == Approx(2.5));
  CHECK(lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  EwaldHandler3D handler(ref, lattice.LR_kc);

  // make sure initBreakup changes the default sigma
  CHECK(handler.Sigma == Approx(lattice.LR_kc));
  handler.initBreakup(ref);
  CHECK(handler.Sigma == Approx(std::sqrt(lattice.LR_kc / (2.0 * lattice.LR_rc))));

  std::cout << "handler.MaxKshell is " << handler.MaxKshell << std::endl;
  CHECK(handler.MaxKshell == 78);
  CHECK(handler.LR_rc == Approx(2.5));
  CHECK(handler.LR_kc == Approx(12));

  mRealType r, dr, rinv;
  mRealType vsr, vlr;
  int nr = 101;
  dr     = 5.0 / nr; // L/[# of grid points]
  for (int ir = 1; ir < nr; ir++)
  {
    r    = ir * dr;
    rinv = 1. / r;
    vsr  = handler.evaluate(r, rinv);
    vlr  = handler.evaluateLR(r);
    // short-range part must vanish after rcut
    if (r > 2.5)
      CHECK(vsr == Approx(0.0));
    // sum must recover the Coulomb potential
    CHECK(vsr + vlr == Approx(rinv));
  }
}

/** evalaute bare Coulomb derivative using EwaldHandler3D
 */
TEST_CASE("ewald3d df", "[lrhandler]")
{
  Lattice lattice;
  lattice.BoxBConds     = true;
  lattice.LR_dim_cutoff = 30.;
  lattice.R.diagonal(5.0);
  lattice.reset();
  CHECK(lattice.Volume == Approx(125));
  lattice.SetLRCutoffs(lattice.Rv);
  //lattice.printCutoffs(app_log());
  CHECK(lattice.LR_rc == Approx(2.5));
  CHECK(lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  EwaldHandler3D handler(ref, lattice.LR_kc);

  // make sure initBreakup changes the default sigma
  CHECK(handler.Sigma == Approx(lattice.LR_kc));
  handler.initBreakup(ref);
  CHECK(handler.Sigma == Approx(std::sqrt(lattice.LR_kc / (2.0 * lattice.LR_rc))));

  std::cout << "handler.MaxKshell is " << handler.MaxKshell << std::endl;
  CHECK(handler.MaxKshell == 78);
  CHECK(handler.LR_rc == Approx(2.5));
  CHECK(handler.LR_kc == Approx(12));

  mRealType r, dr, rinv;
  mRealType rm, rp; // minus (m), plus (p)
  mRealType vsrm, vsrp, dvsr, vlrm, vlrp, dvlr;
  dr                           = 0.00001; // finite difference step
  std::vector<mRealType> rlist = {0.1, 0.5, 1.0, 2.0, 2.5};
  for (auto it = rlist.begin(); it != rlist.end(); ++it)
  {
    r = *it;
    // test short-range piece
    rm   = r - dr;
    rinv = 1. / rm;
    vsrm = handler.evaluate(rm, rinv);
    vlrm = handler.evaluateLR(rm);
    rp   = r + dr;
    rinv = 1. / rp;
    vsrp = handler.evaluate(rp, rinv);
    vlrp = handler.evaluateLR(rp);
    dvsr = (vsrp - vsrm) / (2 * dr);
    rinv = 1. / r;
    CHECK(handler.srDf(r, rinv) == Approx(dvsr));
    // test long-range piece
    dvlr = (vlrp - vlrm) / (2 * dr);
    CHECK(handler.lrDf(r) == Approx(dvlr));
  }
}

} // namespace qmcplusplus
