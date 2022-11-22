#ifndef __EA_EQUATION_OF_MOTION_CCSD_DAVIDSON
#define __EA_EQUATION_OF_MOTION_CCSD_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

class UCcsdEAEquationOfMotionDavidson : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UCcsdEAEquationOfMotionDavidson);
  UCcsdEAEquationOfMotionDavidson(std::vector<Argument> const &argumentList);
  virtual ~UCcsdEAEquationOfMotionDavidson();

  virtual void run();

  template <typename F>
  void run();

protected:
  static constexpr int DEFAULT_MAX_ITERATIONS = 16;
};

} // namespace sisi4s

#endif
