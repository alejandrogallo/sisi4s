#ifndef COMPLEX_TENSOR_CONTRACTION_DEFINED
#define COMPLEX_TENSOR_CONTRACTION_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    ComplexTensorContraction,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new ComplexTensorContraction(argumentList);
    });
} // namespace sisi4s

#endif
