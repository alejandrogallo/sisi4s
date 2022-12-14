#ifndef COMPLEX_TENSOR_WRITER_DEFINED
#define COMPLEX_TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class ComplexTensorWriter : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorWriter);
  ComplexTensorWriter(std::vector<Argument> const &argumentList);
  virtual ~ComplexTensorWriter();
  /**
   * \brief Writes the real tensor data given as Data argument to a file.
   */
  virtual void run();
};
} // namespace sisi4s

#endif
