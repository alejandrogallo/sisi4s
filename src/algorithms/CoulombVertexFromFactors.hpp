#ifndef COULOMB_VERTEX_FROM_FACTORS_DEFINED
#define COULOMB_VERTEX_FROM_FACTORS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/MachineTensor.hpp>
#include <memory>

namespace sisi4s {
/**
 * \brief Caclulates the Coulomb vertex \f$\Gamma^q_{rF}\f$ from the given
 * given factor orbitals \f$\Pi^R_r\f$ and Coulomb factors \f$\Lambda^R_F\f$.
 */
DEFINE_ALGORITHM_HEADER(

    CoulombVertexFromFactors,

    template <typename T, typename MT>
    void run(const bool dryRun););
} // namespace sisi4s

#endif
