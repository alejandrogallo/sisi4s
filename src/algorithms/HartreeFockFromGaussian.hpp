#ifndef _HARTREE_FOCK_FROM_GAUSSIAN_ALGO_DEFINED
#define _HARTREE_FOCK_FROM_GAUSSIAN_ALGO_DEFINED

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace cc4s {

  class HartreeFockFromGaussian: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(HartreeFockFromGaussian);
    HartreeFockFromGaussian(
      std::vector<Argument> const &argumentList
    );
    virtual ~HartreeFockFromGaussian();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };
}

#endif

