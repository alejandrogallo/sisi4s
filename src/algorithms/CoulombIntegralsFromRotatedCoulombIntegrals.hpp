#ifndef _COULOMB_INTEGRALS_FROM_ROTATED_COULOMB_INTEGRALS
#define _COULOMB_INTEGRALS_FROM_ROTATED_COULOMB_INTEGRALS

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

class CoulombIntegralsFromRotatedCoulombIntegrals : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CoulombIntegralsFromRotatedCoulombIntegrals);
  CoulombIntegralsFromRotatedCoulombIntegrals(
      std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~CoulombIntegralsFromRotatedCoulombIntegrals() {}
  virtual void run();
};
} // namespace sisi4s

#endif
