#ifndef Ccsdt_EQUATION_OF_MOTION_DAVIDSON
#define Ccsdt_EQUATION_OF_MOTION_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

  class CcsdtEquationOfMotionDavidson: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdtEquationOfMotionDavidson);
    CcsdtEquationOfMotionDavidson(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdtEquationOfMotionDavidson();

    virtual void run();

    template<typename F>
    void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };

}

#endif

