#ifndef FROM_COMPLEX_TENSOR_DEFINED
#define FROM_COMPLEX_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class FromComplexTensor : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(FromComplexTensor);
  FromComplexTensor(std::vector<Argument> const &argumentList);
  virtual ~FromComplexTensor();
  virtual void run();
};
} // namespace sisi4s

#endif
