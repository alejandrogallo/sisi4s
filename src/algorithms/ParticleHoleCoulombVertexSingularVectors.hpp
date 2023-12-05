#ifndef PARTICLE_HOLE_COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    ParticleHoleCoulombVertexSingularVectors,

    static double constexpr DEFAULT_REDUCTION = 0.5;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;);
} // namespace sisi4s

#endif
