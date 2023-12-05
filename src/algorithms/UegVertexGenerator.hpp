#ifndef UEG_VERTEX_GENERATOR_DEFINED
#define UEG_VERTEX_GENERATOR_DEFINED

#include <array>

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
using ivec = std::array<int, 3>;
using dvec = std::array<double, 4>;

DEFINE_ALGORITHM_HEADER(

    UegVertexGenerator,

    double evalMadelung(double volume);
    double Vijji(const dvec a, const dvec b, const double v);

    bool halfGrid;
    size_t No, Nv, NF;
    double rs, madelung;);
} // namespace sisi4s

#endif
