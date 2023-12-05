#ifndef COMPLEX_TENSOR_SUM_DEFINED
#define COMPLEX_TENSOR_SUM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    ComplexTensorSum,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new ComplexTensorSum(argumentList);
    });
} // namespace sisi4s

#endif
