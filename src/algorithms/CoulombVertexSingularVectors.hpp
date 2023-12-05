#ifndef COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED
#define COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    CoulombVertexSingularVectors,

    static double constexpr DEFAULT_FIELD_VARIABLES_RANK = 0.5;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES_SIZE = -1;);
} // namespace sisi4s

#endif
