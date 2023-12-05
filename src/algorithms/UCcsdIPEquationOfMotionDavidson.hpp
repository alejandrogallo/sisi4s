#ifndef __IP_EQUATION_OF_MOTION_CCSD_DAVIDSON
#define __IP_EQUATION_OF_MOTION_CCSD_DAVIDSON

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    UCcsdIPEquationOfMotionDavidson,

    template <typename F>
    void run();

    static constexpr int DEFAULT_MAX_ITERATIONS = 16;);

} // namespace sisi4s

#endif
