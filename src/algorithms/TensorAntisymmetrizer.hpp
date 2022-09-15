#ifndef ANTISYMMETRIZE_TENSOR_DEFINED
#define ANTISYMMETRIZE_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
  struct TensorAntisymmetrizer: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(TensorAntisymmetrizer);

    TensorAntisymmetrizer(
      std::vector<Argument> const &argumentList): Algorithm(argumentList) {}
    ~TensorAntisymmetrizer(){};

    virtual void run();
  };
}

#endif

