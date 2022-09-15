#ifndef TENSOR_REDUCE_IJIJ_DEFINED
#define TENSOR_REDUCE_IJIJ_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
  struct TensorReduceijij: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(TensorReduceijij);
    TensorReduceijij(
      std::vector<Argument> const &argumentList): Algorithm(argumentList) {}
    ~TensorReduceijij(){};
    virtual void run();
  };
}

#endif

