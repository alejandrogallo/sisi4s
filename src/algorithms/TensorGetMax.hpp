#ifndef TENSOR_GET_MAX_DEFINED 
#define TENSOR_GET_MAX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  struct TensorGetMax: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(TensorGetMax);
    TensorGetMax(
      std::vector<Argument> const &argumentList): Algorithm(argumentList) {}
    ~TensorGetMax(){};
    /**
     * \brief Writes the real tensor data given as Data argument to a file.
     */
    virtual void run();
  };
}

#endif

