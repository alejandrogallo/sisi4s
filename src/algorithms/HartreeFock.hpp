#ifndef _HARTREE_FOCK_ALGO_DEFINED
#define _HARTREE_FOCK_ALGO_DEFINED

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace cc4s {

  class HartreeFock: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(HartreeFock);
    HartreeFock(
      std::vector<Argument> const &argumentList
    );
    virtual ~HartreeFock();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };
}

#endif

