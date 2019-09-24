#ifndef NATURALTRANSITIONORBITALSTRANSFORMATIONFROMPARTICLEHOLEDENSITYMATRIX
#define NATURALTRANSITIONORBITALSTRANSFORMATIONFROMPARTICLEHOLEDENSITYMATRIX

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace cc4s {

  struct NaturalTransitionOrbitalsFromRhoAI: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(NaturalTransitionOrbitalsFromRhoAI);
    NaturalTransitionOrbitalsFromRhoAI(
      std::vector<Argument> const &argumentList
    ): Algorithm(argumentList){};
    virtual ~NaturalTransitionOrbitalsFromRhoAI(){};

    virtual void run();

    template<typename F>
    void run();

  };

}

#endif

