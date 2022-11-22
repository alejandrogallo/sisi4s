#ifndef PSEUDO_INVERSE_SVD_DEFINED
#define PSEUDO_INVERSE_SVD_DEFINED

#include <math/Complex.hpp>
#include <util/Tensor.hpp>
#include <DryTensor.hpp>

#include <util/Tensor.hpp>
#include <random>

namespace sisi4s {
template <typename F>
class PseudoInverseSvd {
public:
  PseudoInverseSvd(Tensor<F> &A, double epsilon = 1e-12);
  CTF::Matrix<F> &get();

protected:
  CTF::Matrix<F> inverse;
};

template <typename F>
class DryPseudoInverseSvd {
public:
  DryPseudoInverseSvd(DryTensor<F> const &matrix);
  DryMatrix<F> &get();

protected:
  DryMatrix<F> inverse;
};
} // namespace sisi4s

#endif
