#ifndef MATH_FUNCTIONS_DEFINED
#define MATH_FUNCTIONS_DEFINED

#include <math/Complex.hpp>
#include <cmath>
#include <string>
#include <util/Tensor.hpp>
#include <util/Log.hpp>

namespace sisi4s {
// constants
template <typename F = double>
constexpr F Pi() {
  return std::acos(F(-1));
}

template <typename F = double>
constexpr F Tau() {
  return 2 * Pi<F>();
}

// univariate functions
template <typename F = double>
inline F sqrt(F const x) {
  return std::sqrt(x);
}

template <typename F = double>
inline F abs(F const x) {
  return std::abs(x);
}

template <typename F>
inline F conj(F const x) {
  return std::conj(x);
}

template <>
inline double conj(double const x) {
  return x;
}

// bivariate functions
template <typename F = double>
inline F dot(F const x, F const y) {
  return conj(x) * y;
}

/**
 * \brief Calculates only the real part of x*conj(y).
 */
template <typename F = double>
inline F realDot(F const x, F const y) {
  return std::real(x * conj(y));
}

template <typename F = double>
inline F divide(F const x, F const y) {
  return x / y;
}

template <typename F = double>
inline F multiply(F const x, F const y) {
  return x * y;
}

template <typename F>
inline double frobeniusNorm(Tensor<F> &t) {
  char *indices(new char[t.order + 1]);
  for (int index(0); index < t.order; ++index) indices[index] = 'a' + index;
  indices[t.order] = 0;
  CTF::Bivar_Function<F> fRealDot(&sisi4s::realDot<F>);
  CTF::Scalar<F> s(*t.wrld);
  s.contract(1.0, t, indices, t, indices, 0.0, "", fRealDot);
  return std::sqrt(std::real(s.get_val()));
}

/**
 * \brief Apply a permutation operator and antisymmetrize accordingly, e.g.
 *  antiSymmetrize(X, "abcd", "abdc") is replaced by
 *  -1 X["abcd"] = X["abcd"] + sign(abcd -> abdc) X["abdc"]
 *
 */
template <typename F>
inline void symmetrize(std::string indices,
                       std::string permuted,
                       Tensor<F> &t,
                       F prefactor = 1) {
  t[indices.c_str()] += prefactor * t[permuted.c_str()];
}
template <typename F>
inline void checkAntisymmetry(Tensor<F> &t) {
  Tensor<F> testResultUp(t);
  Tensor<F> testResultDown(t);
  F normValue;
  testResultUp["abij"] += testResultUp["baij"];
  testResultDown["abij"] += testResultDown["abji"];
  normValue = testResultUp.norm1();
  normValue += testResultDown.norm1();
  if (normValue >= 1e-3) {
    t.print();
    LOG(0, "AntisymmetryCheck")
        << t.get_name() << ": zero tensor norm " << normValue << std::endl;
    exit(1);
  }
}
} // namespace sisi4s

#endif
