/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_NATURAL_ORBITALS
#define MP2_NATURAL_ORBITALS

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates MP2 natural orbitals
 */
class Mp2NaturalOrbitals : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(Mp2NaturalOrbitals);
  Mp2NaturalOrbitals(std::vector<Argument> const &argumentList);
  virtual ~Mp2NaturalOrbitals();
  /**
   * \brief Calculates MP2 energy from Coulomb integrals Vabij
   */
  virtual void run();

protected:
};
} // namespace sisi4s

#endif
