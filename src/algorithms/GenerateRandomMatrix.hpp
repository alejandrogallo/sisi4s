#ifndef GENERATE_RANDOM_MATRIX_DEFINED
#define GENERATE_RANDOM_MATRIX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    GenerateRandomMatrix,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomMatrix(argumentList);
    });
} // namespace sisi4s

#endif
