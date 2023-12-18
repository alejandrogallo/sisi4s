#ifndef PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED
#define PARTICLE_HOLE_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates the Coulomb Integrals \f$V_{ij}^{ab} from the Particle Hole
 * Coulomb Vertex \f$\Gamma_{iG}^a\f$ and stores them in a CTF Tensor Vabij
 * The argument of the integrals is PPHHCoulombIntegrals.
 */
DEFINE_ALGORITHM_HEADER(

    ParticleHoleCoulombIntegrals,

    /** \brief The Coulomb Vertex GammaGai  */
    Tensor<sisi4s::complex> *GammaGai;
    /** \brief The Coulomb integrals Vabij  */
    Tensor<double> * Vabij;);
} // namespace sisi4s

#endif
