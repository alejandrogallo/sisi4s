#ifndef TENSOR_CONTRACTION_DEFINED
#define TENSOR_CONTRACTION_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    TensorContraction,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new TensorContraction(argumentList);
    });
} // namespace sisi4s

#endif
