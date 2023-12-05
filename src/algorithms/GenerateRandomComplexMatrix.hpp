#ifndef GENERATE_RANDOM_COMPLEX_MATRIX_DEFINED
#define GENERATE_RANDOM_COMPLEX_MATRIX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    GenerateRandomComplexMatrix,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomComplexMatrix(argumentList);
    });
} // namespace sisi4s

#endif
