#ifndef QMCPLUSPLUS_LENNARDJONESPBC_H
#define QMCPLUSPLUS_LENNARDJONESPBC_H
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {
    struct LennardJonesPBC: public QMCHamiltonianBase {
        int nParticles;
        DistanceTableData* d_table;
        // SymmetricDTD<RealType, OHMMS_DIM, SUPERCELL_BULK>* d_table;
        ParticleSet &HeConfig;
        Return_t rcut, sigma, epsilon;
        RealType tailCorrection;        //, v2_shift;

        LennardJonesPBC(ParticleSet &atoms, RealType sigma_rhs = 4.830139974, RealType epsilon_rhs = 3.236460251e-5);
        ~LennardJonesPBC() { }

        void resetTargetParticleSet(ParticleSet& P);

        Return_t evaluate(ParticleSet& P);

        /// Frenkel & Smit (1e) Eq.(3.2.5) ... without the rho/2 factor
        Return_t v2_tail(RealType r) {
          RealType srinv3 = std::pow(sigma/r, 3.0);
          return (16.0*M_PI/3.0)*epsilon*sigma*sigma*sigma*((srinv3*srinv3*srinv3/3.0) - srinv3);
        }

        Return_t v2(RealType r) {
          RealType srinv6 = std::pow(sigma/r, 6.0);
          return 4.0*epsilon*(srinv6*srinv6 - srinv6);
        }
        // Return_t v2_rscaled(RealType srinv6) { return 4.0*epsilon*(srinv6*srinv6 - srinv6); }

        // Return_t v2_ts(RealType r) { return v2(r) - v2_shift; } // truncated and shifted

        /** Do nothing */
        bool put(xmlNodePtr cur) {
            return true;
        }

        bool get(std::ostream& os) const {
            os << "LennardJonesPBC: " << HeConfig.getName();
            return true;
        }

        QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi) {
            return new LennardJonesPBC(qp);
        }
    };
}
#endif
