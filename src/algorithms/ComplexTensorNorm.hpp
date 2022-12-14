#ifndef COMPLEX_TENSOR_NORM_DEFINED
#define COMPLEX_TENSOR_NORM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class ComplexTensorNorm : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorNorm);
  ComplexTensorNorm(std::vector<Argument> const &argumentList);
  virtual ~ComplexTensorNorm();
  virtual void run();

  static Algorithm *create(std::vector<Argument> const &argumentList) {
    return new ComplexTensorNorm(argumentList);
  }
};
} // namespace sisi4s

#endif
