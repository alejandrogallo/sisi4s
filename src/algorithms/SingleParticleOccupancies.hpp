#ifndef SINGLE_PARTICLE_OCCUPANCIES_DEFINED
#define SINGLE_PARTICLE_OCCUPANCIES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Evaluates
 * \f$\langle\Psi|\hat N_p|\Psi\rangle/\langle\Psi|\Psi\rangle\f$
 * given the DoublesAmplitudes from a linearized coupled cluster theory.
 */
DEFINE_ALGORITHM_HEADER(

    SingleParticleOccupancies,

);
} // namespace sisi4s

#endif
