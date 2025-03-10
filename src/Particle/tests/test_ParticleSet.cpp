//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
TEST_CASE("ParticleSet distance table management", "[particle]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elecs(simulation_cell);

  ions.setName("ions");
  elecs.setName("electrons");

  const int ii_table_id = ions.addTable(ions);
  const int ie_table_id = ions.addTable(elecs);
  const int ei_table_id = elecs.addTable(ions);
  const int ee_table_id = elecs.addTable(elecs);

  REQUIRE(ii_table_id == 0);
  REQUIRE(ie_table_id == 1);
  REQUIRE(ei_table_id == 0);
  REQUIRE(ee_table_id == 1);

  // second query
  const int ii_table_id2 = ions.addTable(ions);
  const int ie_table_id2 = ions.addTable(elecs);
  const int ei_table_id2 = elecs.addTable(ions);
  const int ee_table_id2 = elecs.addTable(elecs);

  REQUIRE(ii_table_id2 == 0);
  REQUIRE(ie_table_id2 == 1);
  REQUIRE(ei_table_id2 == 0);
  REQUIRE(ee_table_id2 == 1);

  REQUIRE(&(ions.getDistTable(ii_table_id2).get_origin()) == &ions);
  REQUIRE(&(ions.getDistTable(ie_table_id2).get_origin()) == &elecs);
  REQUIRE(&(elecs.getDistTable(ei_table_id2).get_origin()) == &ions);
  REQUIRE(&(elecs.getDistTable(ee_table_id2).get_origin()) == &elecs);

  ParticleSet elecs_copy(elecs);
  REQUIRE(elecs_copy.getDistTable(ei_table_id2).get_origin().getName() == "ions");
  REQUIRE(elecs_copy.getDistTable(ee_table_id2).get_origin().getName() == "electrons");
}

TEST_CASE("symmetric_distance_table OpenBC", "[particle]")
{
  const SimulationCell simulation_cell;
  ParticleSet source(simulation_cell);

  source.setName("electrons");

  source.create({2});
  source.R[0] = {0.0, 1.0, 2.0};
  source.R[1] = {1.1, 1.0, 3.2};
  source.update();
  /// make sure getCoordinates().getAllParticlePos() is updated no matter SoA or AoS.
  CHECK(source.getCoordinates().getAllParticlePos()[0][1] == Approx(1.0));
  CHECK(source.getCoordinates().getAllParticlePos()[1][2] == Approx(3.2));

  const int TableID = source.addTable(source);
  source.update();
  const auto& d_aa      = source.getDistTableAA(TableID);
  const auto& aa_dists  = d_aa.getDistances();
  const auto& aa_displs = d_aa.getDisplacements();

  CHECK(aa_dists[1][0] == Approx(1.62788206));
  CHECK(aa_displs[1][0][0] == Approx(-1.1));
  CHECK(aa_displs[1][0][1] == Approx(0.0));
  CHECK(aa_displs[1][0][2] == Approx(-1.2));
}

TEST_CASE("symmetric_distance_table PBC", "[particle]")
{
  Lattice lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R = ParticleSet::Tensor_t(6.74632230, 6.74632230, 0.00000000, 0.00000000, 3.37316115, 3.37316115, 3.37316115,
                                    0.00000000, 3.37316115);
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet source(simulation_cell);

  source.setName("electrons");

  source.create({4});
  source.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
  source.R[1] = ParticleSet::PosType(1.68658058, 1.68658058, 1.68658058);
  source.R[2] = ParticleSet::PosType(3.37316115, 3.37316115, 0.00000000);
  source.R[3] = ParticleSet::PosType(5.05974172, 5.05974172, 1.68658058);

  const int TableID = source.addTable(source);
  source.update();
  const auto& d_aa      = source.getDistTableAA(TableID);
  const auto& aa_dists  = d_aa.getDistances();
  const auto& aa_displs = d_aa.getDisplacements();

  CHECK(aa_dists[2][1] == Approx(2.9212432441));
  CHECK(aa_displs[2][1][0] == Approx(-1.68658057));
  CHECK(aa_displs[2][1][1] == Approx(-1.68658057));
  CHECK(aa_displs[2][1][2] == Approx(1.68658057));
}

TEST_CASE("particle set lattice with vacuum", "[particle]")
{
  // PPP case
  Lattice lattice;
  lattice.BoxBConds          = true;
  lattice.R                  = {1.0, 2.0, 3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  lattice.explicitly_defined = true;

  lattice.VacuumScale = 2.0;
  lattice.reset();
  {
    const SimulationCell simulation_cell(lattice);
    ParticleSet source(simulation_cell);
    source.setName("electrons");
    source.createSK();

    CHECK(lattice.SuperCellEnum == SUPERCELL_BULK);
    CHECK(source.getLRBox().R(0, 0) == 1.0);
    CHECK(source.getLRBox().R(0, 1) == 2.0);
    CHECK(source.getLRBox().R(0, 2) == 3.0);
  }

  // PPN case
  lattice.BoxBConds[2] = false;
  lattice.reset();
  {
    const SimulationCell simulation_cell(lattice);
    ParticleSet source(simulation_cell);
    source.setName("electrons");
    source.createSK();

    CHECK(lattice.SuperCellEnum == SUPERCELL_SLAB);
    CHECK(source.getLRBox().R(2, 0) == 0.0);
    CHECK(source.getLRBox().R(2, 1) == 0.0);
    CHECK(source.getLRBox().R(2, 2) == 2.0);
  }

  // PNN case
  lattice.BoxBConds[1] = false;
  lattice.reset();
  {
    const SimulationCell simulation_cell(lattice);
    ParticleSet source(simulation_cell);
    source.setName("electrons");
    source.createSK();

    CHECK(lattice.SuperCellEnum == SUPERCELL_WIRE);
    CHECK(source.getLRBox().R(0, 0) == 1.0);
    CHECK(source.getLRBox().R(0, 1) == 2.0);
    CHECK(source.getLRBox().R(0, 2) == 3.0);
    CHECK(source.getLRBox().R(1, 0) == 0.0);
    CHECK(source.getLRBox().R(1, 1) == 2.0);
    CHECK(source.getLRBox().R(1, 2) == 0.0);
    CHECK(source.getLRBox().R(2, 0) == 0.0);
    CHECK(source.getLRBox().R(2, 1) == 0.0);
    CHECK(source.getLRBox().R(2, 2) == 2.0);

    source.create({1});
    source.R[0] = {0.5, 1.5, 2.0};

    ParticleSet source_clone(source);
    RefVectorWithLeader<ParticleSet> p_ref_list(source, {source, source_clone});

    ResourceCollection pset_res("test_pset_res");
    source.createResource(pset_res);
    ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_ref_list);

    std::vector<ParticleSet::PosType> displs{{1.0, 2.0, 3.0}, {0.0, 1.0, 0.0}};
    std::vector<bool> are_valid(displs.size());
    ParticleSet::mw_makeMove(p_ref_list, 0, displs, are_valid);
    CHECK(are_valid[0]);  // "p" direction, outside box allowed
    CHECK(!are_valid[1]); // "n" direction, outside box not allowed
  }
}

} // namespace qmcplusplus
