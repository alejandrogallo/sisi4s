#ifndef COMPLEX_TENSOR_NORM_DEFINED
#define COMPLEX_TENSOR_NORM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    ComplexTensorNorm,

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new ComplexTensorNorm(argumentList);
    });

} // namespace sisi4s

#endif
