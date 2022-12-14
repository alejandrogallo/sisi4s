#ifndef PSEUDO_INVERSE_HERMITIAN_SVD_DEFINED
#define PSEUDO_INVERSE_HERMITIAN_SVD_DEFINED

#include <math/Complex.hpp>
#include <DryTensor.hpp>

#include <util/Tensor.hpp>
#include <random>

namespace sisi4s {
template <typename F>
class PseudoInverseHermitianSvd {
public:
  PseudoInverseHermitianSvd(CTF::Matrix<F> &matrix, double epsilon = 1 - 12);
  CTF::Matrix<F> &get();

protected:
  CTF::Matrix<F> inverse;
};

template <typename F>
class DryPseudoInverseHermitianSvd {
public:
  DryPseudoInverseHermitianSvd(DryMatrix<F> const &matrix);
  DryMatrix<F> &get();

protected:
  DryMatrix<F> inverse;
};
} // namespace sisi4s

#endif
