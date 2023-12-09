#ifndef UEG_VERTEX_GENERATOR_DEFINED
#define UEG_VERTEX_GENERATOR_DEFINED

#include <array>

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
using ivec = std::array<int, 3>;
using dvec = std::array<double, 4>;

DEFINE_ALGORITHM_HEADER(

    UegVertexGenerator,

    bool halfGrid;

);
} // namespace sisi4s

#endif
