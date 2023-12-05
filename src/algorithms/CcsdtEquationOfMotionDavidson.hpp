#ifndef Ccsdt_EQUATION_OF_MOTION_DAVIDSON
#define Ccsdt_EQUATION_OF_MOTION_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    CcsdtEquationOfMotionDavidson,

    template <typename F>
    void run();

    static constexpr int DEFAULT_MAX_ITERATIONS = 16;);

} // namespace sisi4s

#endif
