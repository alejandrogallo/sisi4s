#ifndef CCSD_DIAGRAMMATIC_DECOMPOSITION
#define CCSD_DIAGRAMMATIC_DECOMPOSITION

#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <vector>
#include <string>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    CcsdDiagrammaticDecomposition,

    void evaluateEnergy(std::string diagramType,
                        Tensor<double> &deltaabij,
                        Tensor<double> &deltaai,
                        Tensor<double> &Rabij,
                        Tensor<double> &Rai);
    void sliceIntoResiduum(Tensor<double> &Rxyij,
                           int a,
                           int b,
                           Tensor<double> &Rabij););
} // namespace sisi4s
#endif
