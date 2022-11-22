#ifndef __CCSD_EQUATION_OF_MOTION_DAVIDSON
#define __CCSD_EQUATION_OF_MOTION_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

class CcsdEquationOfMotionDavidson : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdEquationOfMotionDavidson);
  CcsdEquationOfMotionDavidson(std::vector<Argument> const &argumentList);
  virtual ~CcsdEquationOfMotionDavidson();

  virtual void run();

  template <typename F>
  void run();

protected:
  static constexpr int DEFAULT_MAX_ITERATIONS = 16;
};

} // namespace sisi4s

#endif
