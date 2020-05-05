#ifndef TENOSR_SLICER_GOD_IS_IN_THE_DOLLARS
#define TENOSR_SLICER_GOD_IS_IN_THE_DOLLARS

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  struct TensorSlicer: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(TensorSlicer);
    TensorSlicer(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
    ~TensorSlicer(){};
    virtual void run();
  };
}

#endif

