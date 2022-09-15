#ifndef UNRESTRICT_TENSOR_DEFINED
#define UNRESTRICT_TENSOR_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
  struct TensorUnrestricter: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(TensorUnrestricter);

    TensorUnrestricter(
      std::vector<Argument> const &argumentList): Algorithm(argumentList) {}
    ~TensorUnrestricter(){};

    virtual void run();
  };
}

#endif

