#ifndef TENSOR_GET_MAX_DEFINED
#define TENSOR_GET_MAX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class TensorGetMax : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(TensorGetMax);
  TensorGetMax(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~TensorGetMax() {}
  /**
   * \brief Writes the real tensor data given as Data argument to a file.
   */
  virtual void run();
};
} // namespace sisi4s

#endif
