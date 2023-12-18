#ifndef FINITE_SIZE_CORRECTION_DEFINED
#define FINITE_SIZE_CORRECTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Interpolation.hpp>
#include <math/Vector.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    FiniteSizeCorrection,

    int NG;
    std::vector<double> GLengths;
    std::vector<double> averageSGs;
    std::vector<double> meanErrorSG;
    std::vector<double> structureFactors;
    std::vector<double> VofG;
    double GC;
    double inter3D;
    double sum3D;
    class Momentum;
    Momentum * fibonacciGrid;
    Momentum * cartesianGrid;
    void readFromFile();
    void calculateRealStructureFactor();
    //    template <typename F>
    //    F calculateStructureFactor(Tensor<F> &Vabij);

    void calculateComplexStructureFactor();

    void constructFibonacciGrid(double R, int N);
    void interpolation3D();
    bool IsInSmallBZ(Vector<double> point,
                     double scale,
                     std::vector<sisi4s::Vector<double>> smallBZ);
    double SGxVG(sisi4s::Inter1D<double> Int1d, double x);
    double integrate(sisi4s::Inter1D<double> Int1d,
                     double start,
                     double end,
                     int steps);
    double simpson(sisi4s::Inter1D<double> Int1d, double x, double h);
    void calculateFiniteSizeCorrection();

    void dryCalculateStructureFactor();
    void dryInterpolation3D();
    void dryCalculateFiniteSizeCorrection();

    void extrapolation(double minG, double maxG, int basisSetExtrapolation);
    double simplestWindow(double Gmin, double Gmax, double G);
    double integrateSimplestWindow(double Gmin, double Gmax);
    double leastSquareFit(std::vector<double> fitabsG,
                          std::vector<double> fitSF);

    void basisSetCompleteness(););

} // namespace sisi4s

#endif
