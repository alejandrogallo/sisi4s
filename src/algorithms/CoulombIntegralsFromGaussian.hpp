#ifndef _COULOMB_INTEGRALS_FROM_GAUSSIAN
#define _COULOMB_INTEGRALS_FROM_GAUSSIAN

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

class CoulombIntegralsFromGaussian : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CoulombIntegralsFromGaussian);
  CoulombIntegralsFromGaussian(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~CoulombIntegralsFromGaussian() {}
  virtual void run();
};
} // namespace sisi4s

#endif
