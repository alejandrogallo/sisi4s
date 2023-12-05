#ifndef BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED
#define BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED

#include <vector>
#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    BasisSetExtrapolation,

    void evaluateQGG(int orbitalPairStart, int orbtialPairEnd, int slice);
    void fitF12(int type, real minG, real maxG);
    void calculateNewSF(int type,
                        real gamma,
                        Tensor<double> *coulombKernel,
                        Tensor<double> *newSF,
                        Tensor<double> *resNewSF);
    void invertQGG(););
} // namespace sisi4s

#endif
