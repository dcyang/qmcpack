#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCHamiltonians/LennardJonesPBC.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
    LennardJonesPBC::LennardJonesPBC(ParticleSet &He_atoms, RealType sigma_rhs, RealType epsilon_rhs):
        HeConfig(He_atoms), d_table(0), rcut(He_atoms.LRBox.LR_rc), sigma(sigma_rhs), epsilon(epsilon_rhs)
    {
        //set the distance tables
        d_table = DistanceTable::add(He_atoms,DT_AOS);  // AOS = array of structs
        // d_table = SymmetricDTD<RealType, OHMMS_DIM, SUPERCELL_BULK>::add(He_atoms,DT_AOS);  // AOS = array of structs
        nParticles = He_atoms.getTotalNum();
        // v2_shift = v2(rcut);
        tailCorrection = (nParticles*nParticles/(16.0*rcut))*v2_tail(rcut);
        // coefficient works out to N*rho/2 ... since V = (2*rcut)^3

        set_energy_domain(potential);
        set_quantum_domain(quantum);
    }

    void LennardJonesPBC::resetTargetParticleSet(ParticleSet& P)
    {
        d_table = DistanceTable::add(HeConfig,P,DT_AOS);
        // d_table = SymmetricDTD<RealType, OHMMS_DIM, SUPERCELL_BULK>::add(HeConfig,P,DT_AOS);
        rcut = P.LRBox.LR_rc;
    }

    LennardJonesPBC::Return_t LennardJonesPBC::evaluate(ParticleSet& P)
    {
      //calculate the Electron-Core Dipole matrix
      //CoreElDipole=0.0;
      RealType e = 0.0;
      for (int i=0; i<nParticles; ++i) {
        // XXX: d_table == P.DistTables[0] ... no need to define as separate d_aa
        // XXX: ... need a closer look at M, J, IJ, PairID in d_table
        //      take CoulombPBCAA::evaluate() for example
        for (int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn) {
          if (d_table->r(nn) > rcut) continue;
          // looks like "sphere" carving isn't fully implemented (supporting only Coulomb?)
          // -> doing this manually

          // int j(d_table->J[nn]);
          RealType srinv6 = std::pow(sigma*d_table->rinv(nn), 6.0);
          e += srinv6*srinv6 - srinv6;
          /*
#pragma omp master
          {
            app_log() << "d_table->M[" << i << "] (" << nn << ", " << d_table->J[nn] << ") = " << d_table->r(nn) << "; " << rcut << std::endl;
          }
          */
        }
      }
      e *= 4.0*epsilon;

      e += tailCorrection;

      return Value=e;
    }

    void HamiltonianFactory::addAtomicDimerPotential(xmlNodePtr cur) {
        typedef QMCHamiltonian::Return_t Return_t;
        std::string targetInp(targetPtcl->getName());
        std::string sourceInp(targetPtcl->getName());
        std::string title("AtomAtom");  // ,pbc("yes");
        // std::string forces("no");    // never heard of this being done anywhere
        bool physical = true;
        bool doForce = false;

        double sigma = 0.0, epsilon = 0.0;      // as defined in McMillan, PR 138, A442 (1965)

        OhmmsAttributeSet hAttrib;
        hAttrib.add(title,"id");
        hAttrib.add(title,"name");
        hAttrib.add(targetInp,"target");
        hAttrib.add(sourceInp,"source");
        // hAttrib.add(pbc,"pbc"); // assumed true (yes)
        hAttrib.add(physical,"physical");
        // hAttrib.add(forces,"forces");   

        hAttrib.add(sigma,"sigma");
        hAttrib.add(epsilon,"epsilon");

        hAttrib.put(cur);

        // bool applyPBC = (PBCType && pbc=="yes");
        // bool doForces = (forces == "yes") || (forces == "true");
        // ParticleSet *ptclA=targetPtcl;

        // set to McMillan values if sigma, epsilon are non-positive
        if (sigma < 1.0e-10) sigma = 4.830139974; // bohr radii
        if (epsilon < 1.0e-10) epsilon = 3.236460251e-5; // Ha

        /*
        // debug lines
        app_log() << "<NENE>" << std::endl;
        app_log() << "title = " << title << std::endl;
        app_log() << "targetInp = " << targetInp << std::endl;
        app_log() << "sourceInp = " << sourceInp << std::endl;
        // app_log() << "pbc = " << pbc << std::endl;
        app_log() << "sigma = " << sigma << std::endl;
        app_log() << "epsilon = " << epsilon << std::endl;
        // app_log() << "physical = " << physical << std::endl;
        app_log() << "</NENE>" << std::endl;
        app_log().flush();
        */

        LennardJonesPBC *ljp = new LennardJonesPBC(*targetPtcl, sigma, epsilon);
        ljp->sigma = sigma;
        ljp->epsilon = epsilon;

        targetH->addOperator(ljp, title);
        // appends to the vector targetH->H
    }
}
