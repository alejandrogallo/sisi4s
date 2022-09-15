#ifndef __IP_EQUATION_OF_MOTION_CCSD_DAVIDSON
#define __IP_EQUATION_OF_MOTION_CCSD_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

  class UCcsdIPEquationOfMotionDavidson: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UCcsdIPEquationOfMotionDavidson);
    UCcsdIPEquationOfMotionDavidson(std::vector<Argument> const &argumentList);
    virtual ~UCcsdIPEquationOfMotionDavidson();

    virtual void run();

    template<typename F>
    void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };

}

#endif

