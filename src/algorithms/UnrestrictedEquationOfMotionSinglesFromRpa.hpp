#ifndef __RPA_SINGLES_EQUATION_OF_MOTION_DAVIDSON
#define __RPA_SINGLES_EQUATION_OF_MOTION_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

class UnrestrictedEquationOfMotionSinglesFromRpa : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UnrestrictedEquationOfMotionSinglesFromRpa);
  UnrestrictedEquationOfMotionSinglesFromRpa(
      std::vector<Argument> const &argumentList);
  virtual ~UnrestrictedEquationOfMotionSinglesFromRpa();

  virtual void run();

  template <typename F>
  void run();

protected:
  static constexpr int DEFAULT_MAX_ITERATIONS = 16;
};

} // namespace sisi4s

#endif
