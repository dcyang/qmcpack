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

        LennardJonesPBC(ParticleSet &atoms);
        ~LennardJonesPBC() { }

        void resetTargetParticleSet(ParticleSet& P);

        Return_t evaluate(ParticleSet& P);

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
