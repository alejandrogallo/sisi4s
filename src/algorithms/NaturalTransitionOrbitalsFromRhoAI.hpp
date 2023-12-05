#ifndef NATURALTRANSITIONORBITALSTRANSFORMATIONFROMPARTICLEHOLEDENSITYMATRIX
#define NATURALTRANSITIONORBITALSTRANSFORMATIONFROMPARTICLEHOLEDENSITYMATRIX

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    NaturalTransitionOrbitalsFromRhoAI,

    template <typename F>
    void buildTransformations(Tensor<F> &rho, const std::string name);

    template <typename F>
    void run(););

} // namespace sisi4s

#endif
