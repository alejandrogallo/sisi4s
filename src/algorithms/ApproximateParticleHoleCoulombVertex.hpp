#ifndef APPROXIMATE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED
#define APPROXIMATE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
/* @WIKI
 * \brief Approximates the particle-hole Coulomb vertex $\tilde\Gamma^a_{iG}$
 * using the given set of left singular vectors $U^F_G$ associated to the
 * $\Gamma^a_{iF} = {U^\ast}^G_F \tilde\Gamma^a_{iG}$.
 * largest singular values, The approximated particle-hole Coulomb vertex is
 * given by \param[in] ParticleHoleCoulombVertex (complex tensor) (none) The
 * particle-hole Coulomb vertex $\tilde\Gamma$ to approximate. \param[in]
 * ParticleHoleCoulomgVertexSingularVectors tensor none The left singular
 * vectors $U$ to be used in the transformation. \param[out]
 * ApproximatedParticleHoleCoulombVertex tensor The approximated particle-hole
 * Coulomb vertex $\Gamma$.
 */
DEFINE_ALGORITHM_HEADER(ApproximateParticleHoleCoulombVertex, );
} // namespace sisi4s

#endif
