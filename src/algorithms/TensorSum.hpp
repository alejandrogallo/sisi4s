#ifndef TENSOR_SUM_DEFINED
#define TENSOR_SUM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class TensorSum : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(TensorSum);
  TensorSum(std::vector<Argument> const &argumentList);
  virtual ~TensorSum();
  virtual void run();

  static Algorithm *create(std::vector<Argument> const &argumentList) {
    return new TensorSum(argumentList);
  }
};
} // namespace sisi4s

#endif
