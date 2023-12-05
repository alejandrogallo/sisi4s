#ifndef APPROXIMATE_COULOMB_VERTEX_DEFINED
#define APPROXIMATE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
/* @WIKI
 * \brief Approximates the Coulomb vertex $\tilde\Gamma^q_{rG}$
 * using the given set of left singular vectors $U^F_G$ associated to the
 * $\Gamma^q_{rF} = {U^\ast}^G_F \tilde\Gamma^q_{rG}$.
 * largest singular values, The approximated Coulomb vertex is given by
 * \param[in] CoulombVertex (complex tensor) (none)
 *   The Coulomb vertex $\tilde\Gamma$ to approximate.
 * \param[in] EnergyMatrixTransform tensor none
 *   The left singular vectors $U$ to be used in the transformation.
 * \param[out] ApproximatedCoulombVertex tensor
 *   The approximated Coulomb vertex $\Gamma$.
 */
DEFINE_ALGORITHM_HEADER(ApproximateCoulombVertex, );
} // namespace sisi4s

#endif
