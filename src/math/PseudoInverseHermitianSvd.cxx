#include <math/PseudoInverseHermitianSvd.hpp>

#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>
#include <memory>
#include <util/Tensor.hpp>

using namespace sisi4s;
using std::make_shared;

template <typename F>
PseudoInverseHermitianSvd<F>::PseudoInverseHermitianSvd(CTF::Matrix<F> &A,
                                                        double epsilon)
    : inverse(A) {
  // convert CTF matrices into ScaLapack matrices
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  auto scaA(make_shared<ScaLapackMatrix<F>>(A, &blacsWorld));
  auto scaU(make_shared<ScaLapackMatrix<F>>(*scaA));

  // solve hermitian eigenvectors problem using ScaLapack
  ScaLapackHermitianEigenSystemDc<F> eigenSystem(scaA, scaU);
  double *lambda(new double[scaA->lens[0]]);
  eigenSystem.solve(lambda);

  CTF::Vector<F> D(A.lens[0], *A.wrld, "Lambda");
  int localLambdaCount(A.wrld->rank == 0 ? A.lens[0] : 0);
  F *lambdaValues(new F[localLambdaCount]);
  int64_t *lambdaIndices(new int64_t[localLambdaCount]);
  bool smallLambda(false);
  for (int64_t i(0); i < localLambdaCount; ++i) {
    lambdaIndices[i] = i;
    smallLambda |= std::abs(lambda[i]) < 1e3 * epsilon;
    lambdaValues[i] = (std::abs(lambda[i]) > epsilon) ? 1 / lambda[i] : 0;
  }
  if (smallLambda) {
    LOG(1, "PseudoInverseHermitianSvd")
        << "WARNING: Very small eigenvalues occurred" << std::endl;
  }
  D.write(localLambdaCount, lambdaIndices, lambdaValues);
  CTF::Matrix<F> U(A);
  scaU->write(U);
  auto conjU(U);
  conjugate(conjU);
  delete[] lambdaIndices;
  delete[] lambdaValues;
  delete[] lambda;

  // compose the pseudo inverse from S and the unitary transforms
  inverse["ij"] = U["ik"] * D["k"] * conjU["jk"];
}

template <typename F>
CTF::Matrix<F> &PseudoInverseHermitianSvd<F>::get() {
  return inverse;
}

// instantiate
template PseudoInverseHermitianSvd<sisi4s::Float64>::PseudoInverseHermitianSvd(
    CTF::Matrix<sisi4s::Float64> &matrix,
    double epsilon);
template CTF::Matrix<sisi4s::Float64> &
PseudoInverseHermitianSvd<sisi4s::Float64>::get();

template PseudoInverseHermitianSvd<sisi4s::Complex64>::
    PseudoInverseHermitianSvd(CTF::Matrix<sisi4s::Complex64> &matrix,
                              double epsilon);
template CTF::Matrix<sisi4s::Complex64> &
PseudoInverseHermitianSvd<sisi4s::Complex64>::get();

template <typename F>
DryPseudoInverseHermitianSvd<F>::DryPseudoInverseHermitianSvd(
    DryMatrix<F> const &matrix_)
    : inverse(matrix_) {}

template <typename F>
DryMatrix<F> &DryPseudoInverseHermitianSvd<F>::get() {
  return inverse;
}

// instantiate
template DryPseudoInverseHermitianSvd<sisi4s::Float64>::
    DryPseudoInverseHermitianSvd(DryMatrix<sisi4s::Float64> const &matrix);
template DryMatrix<sisi4s::Float64> &
DryPseudoInverseHermitianSvd<sisi4s::Float64>::get();

template DryPseudoInverseHermitianSvd<sisi4s::Complex64>::
    DryPseudoInverseHermitianSvd(DryMatrix<sisi4s::Complex64> const &matrix);
template DryMatrix<sisi4s::Complex64> &
DryPseudoInverseHermitianSvd<sisi4s::Complex64>::get();
