#ifndef PARTICLE_HOLE_COULOMB_VERTEX_FROM_FACTORS_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_FROM_FACTORS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/MachineTensor.hpp>
#include <memory>

namespace sisi4s {
/**
 * \brief Caclulates the particle hole Coulomb vertex \f$\Gamma^a_{iF}\f$
 * from the given given particle factors orbitals \f$\Pi^R_a\f$,
 * hole factor orbitals \f$\Pi^R_i\f$ and Coulomb factors \f$\Lambda^R_F\f$.
 */
DEFINE_ALGORITHM_HEADER(

    ParticleHoleCoulombVertexFromFactors,

    template <typename T, typename MT>
    void run(const bool dryRun););
} // namespace sisi4s

#endif
