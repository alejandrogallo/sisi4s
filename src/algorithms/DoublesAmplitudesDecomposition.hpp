#ifndef DOUBLES_AMPLITUDES_DECOMPOSITION_DEFINED
#define DOUBLES_AMPLITUDES_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <memory>

namespace sisi4s {
DEFINE_ALGORITHM_HEADER(

    DoublesAmplitudesDecomposition,

    static double constexpr DEFAULT_REDUCTION = 2.0;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;

    void diagonlizeAmplitudes();
    void sliceLargestEigenValues();

    int Nv, No, NvNo, NF, lower, upper;
    double *lambdas;
    complex * sqrtLambdas;
    int64_t lambdasCount, *lambdaIndices;

    std::shared_ptr<Tensor<double>> Taibj, UaiF;
    std::shared_ptr<Tensor<complex>> LFai, sqrtLambdaF;
    Tensor<double> * LambdaF;);
} // namespace sisi4s

#endif
