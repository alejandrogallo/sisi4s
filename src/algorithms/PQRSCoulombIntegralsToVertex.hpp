/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PQRS_COULOMB_INTEGRALS_TO_VERTEX
#define PQRS_COULOMB_INTEGRALS_TO_VERTEX

#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <vector>

namespace sisi4s {
class PQRSCoulombIntegralsToVertex : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(PQRSCoulombIntegralsToVertex);
  PQRSCoulombIntegralsToVertex(std::vector<Argument> const &argumentList);
  virtual ~PQRSCoulombIntegralsToVertex();
  virtual void run();

protected:
};
} // namespace sisi4s

#endif
