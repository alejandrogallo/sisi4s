#ifndef GENERATE_RANDOM_TENSOR_DEFINED
#define GENERATE_RANDOM_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    GenerateRandomTensor,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new GenerateRandomTensor(argumentList);
    });
} // namespace sisi4s

#endif
