#include <array>

#include <math/PseudoInverseSvd.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackSingularValueDecomposition.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

template <typename F>
PseudoInverseSvd<F>::PseudoInverseSvd(Tensor<F> &A, double epsilon)
    : inverse(A.lens[0], A.lens[1]) {
  // convert CTF matrices into ScaLapack matrices
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  // TODO: only works for quadratic matrices
  ScaLapackMatrix<F> ScaA(
      A,
      std::array<int, 2>({(int)A.lens[0], (int)A.lens[1]}).data(),
      &blacsWorld);
  ScaLapackMatrix<F> ScaU(ScaA);
  ScaLapackMatrix<F> ScaVT(ScaA);

  // do SVD using ScaLapack
  ScaLapackSingularValueDecomposition<F> svd(&ScaA, &ScaU, &ScaVT);
  double *sigma(new double[A.lens[0]]);
  svd.decompose(sigma);

  // convert real sigma into complex CTF vector of their pseudo inverses
  // TODO: only works for quadratic matrices
  CTF::Vector<F> S(A.lens[0], *A.wrld, "Sigma");
  int localSigmaCount(A.wrld->rank == 0 ? ScaA.lens[0] : 0);
  F *sigmaValues(new F[localSigmaCount]);
  int64_t *sigmaIndices(new int64_t[localSigmaCount]);
  for (int64_t i(0); i < localSigmaCount; ++i) {
    sigmaIndices[i] = i;
    // invert singular values
    sigmaValues[i] = (sigma[i] > epsilon) ? 1 / sigma[i] : 0;
    // TODO: warn about very small singular values above epsilon
  }
  S.write(localSigmaCount, sigmaIndices, sigmaValues);

  // convert ScaLapack result matrices to CTF
  CTF::Matrix<F> U(A.lens[0], A.lens[1]);
  ScaU.write(U);
  CTF::Matrix<F> VT(A.lens[0], A.lens[1]);
  ScaVT.write(VT);
  delete[] sigmaIndices;
  delete[] sigmaValues;
  delete[] sigma;

  // recompose in CTF to get pseudo inverse matrix
  inverse["ij"] = VT["ki"] * S["k"] * U["jk"];
  CTF::Univar_Function<F> fConj(conj<F>);
  inverse.sum(1.0, inverse, "ij", 0.0, "ij", fConj);
}

template <typename F>
CTF::Matrix<F> &PseudoInverseSvd<F>::get() {
  return inverse;
}

// instantiate
template PseudoInverseSvd<sisi4s::Float64>::PseudoInverseSvd(
    Tensor<sisi4s::Float64> &matrix,
    double epsilon);
template CTF::Matrix<sisi4s::Float64> &PseudoInverseSvd<sisi4s::Float64>::get();

template PseudoInverseSvd<sisi4s::Complex64>::PseudoInverseSvd(
    Tensor<sisi4s::Complex64> &matrix,
    double epsilon);
template CTF::Matrix<sisi4s::Complex64> &
PseudoInverseSvd<sisi4s::Complex64>::get();

template <typename F>
DryPseudoInverseSvd<F>::DryPseudoInverseSvd(DryTensor<F> const &matrix_)
    : inverse(matrix_.lens[0], matrix_.lens[1], NS) {}

template <typename F>
DryMatrix<F> &DryPseudoInverseSvd<F>::get() {
  return inverse;
}

// instantiate
template sisi4s::DryPseudoInverseSvd<sisi4s::Float64>::DryPseudoInverseSvd(
    sisi4s::DryTensor<sisi4s::Float64> const &matrix);

template sisi4s::DryMatrix<sisi4s::Float64> &
sisi4s::DryPseudoInverseSvd<sisi4s::Float64>::get();

template sisi4s::DryPseudoInverseSvd<sisi4s::Complex64>::DryPseudoInverseSvd(
    DryTensor<sisi4s::Complex64> const &matrix);
template sisi4s::DryMatrix<sisi4s::Complex64> &
sisi4s::DryPseudoInverseSvd<sisi4s::Complex64>::get();
