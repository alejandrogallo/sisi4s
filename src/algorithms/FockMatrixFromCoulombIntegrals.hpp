#ifndef FOCK_MATRIX_FROM_COULOMB_INTEGRALS_DEFINED
#define FOCK_MATRIX_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
  struct FockMatrixFromCoulombIntegrals: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(FockMatrixFromCoulombIntegrals);

    FockMatrixFromCoulombIntegrals(
      std::vector<Argument> const &argumentList): Algorithm(argumentList) {}
    ~FockMatrixFromCoulombIntegrals(){};

    virtual void run();
  };
}

#endif

