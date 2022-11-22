#ifndef TENSOR_CONTRACTION_DEFINED
#define TENSOR_CONTRACTION_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class TensorContraction : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(TensorContraction);
  TensorContraction(std::vector<Argument> const &argumentList);
  virtual ~TensorContraction();
  virtual void run();

  static Algorithm *create(std::vector<Argument> const &argumentList) {
    return new TensorContraction(argumentList);
  }
};
} // namespace sisi4s

#endif
