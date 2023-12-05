#include <algorithms/NaturalTransitionOrbitalsFromRhoAI.hpp>
#include <util/SharedPointer.hpp>
#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <util/Log.hpp>
#include <vector>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <numeric> // std::iota
#include <algorithm>

#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(NaturalTransitionOrbitalsFromRhoAI);

IMPLEMENT_EMPTY_DRYRUN(NaturalTransitionOrbitalsFromRhoAI) {}

void NaturalTransitionOrbitalsFromRhoAI::run() {
  if (getIntegerArgument("complex", 1) == 1) {
    run<sisi4s::complex>();
  } else {
    throw new EXCEPTION("No real version of the algo available");
  }
}

template <typename F>
inline void logVectorNorms(LapackMatrix<F> &A, const std::string name) {
  F norm;
  for (int j(0); j < A.getColumns(); j++) {
    norm = F(0);
    for (int i(0); i < A.getRows(); i++) {
      norm += std::conj(A(i, j)) * A(i, j);
    }
    LOG(1, "NaturalTransitionOrbitalsFromRhoAI")
        << "|" << name << "(" << j << ")|**2 " << norm << std::endl;
  }
}

template <typename F>
inline void logOverlap(LapackMatrix<F> &A, const std::string name) {
  F o;
  for (int j(0); j < A.getColumns(); j++) {
    for (int i(0); i < A.getColumns(); i++) {
      o = F(0);
      for (int k(0); k < A.getRows(); k++) {
        o += std::conj(A(k, i)) * A(k, j);
      }
      LOG(1, "NaturalTransitionOrbitalsFromRhoAI")
          << "Overlap_" << name << "(" << i << "," << j << ") " << o
          << std::endl;
    }
  }
}

template <typename F>
void NaturalTransitionOrbitalsFromRhoAI::buildTransformations(
    Tensor<F> &rho,
    const std::string name) {
  // We build first a lapack matrix in order to do the diagonalization
  //
  LapackMatrix<F> rhoMatrix(rho);
  int n(rho.lens[0]);
  int nn[] = {n, n};
  std::vector<int64_t> indices;

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI")
      << "Diagonalizing " << name << " " << n << "x" << n << " block"
      << std::endl;
  LapackGeneralEigenSystem<F> solver(rhoMatrix);

  // get the right eigenvectors
  LapackMatrix<complex> rightEigenVectors(solver.getRightEigenVectors());
  Tensor<F> *rEigenVecs = new Tensor<F>(2, nn, rho.sym, *Sisi4s::world, "r");
  Tensor<F> *overlapMatrix = new Tensor<F>(*rEigenVecs);

  // the indices will hold the vector indices
  if (rEigenVecs->wrld->rank == 0) {
    indices.resize(n * n);
  } else {
    indices.resize(0);
  }
  std::iota(indices.begin(), indices.end(), 0);
  // write the data in the rightEigenVectors into the tensor
  rEigenVecs->write(indices.size(),
                    indices.data(),
                    rightEigenVectors.getValues());
  allocatedTensorArgument<F>(name + "TransformationMatrix", rEigenVecs);

  // log vector norms and overlaps
  logVectorNorms<F>(rightEigenVectors, name);
  // logOverlap<F>(rightEigenVectors, name);

  auto eigenConj = NEW(Tensor<F>, rEigenVecs);
  conjugate(*eigenConj);

  (*overlapMatrix)["ab"] = (*rEigenVecs)["be"] * (*eigenConj)["ea"];
  allocatedTensorArgument<F>(name + "OverlapMatrix", overlapMatrix);

  std::vector<complex> lambdas(solver.getEigenValues());
  Tensor<F> *lambdasTensor =
      new Tensor<F>(1, nn, rho.sym, *Sisi4s::world, "lambdas");
  indices.resize(n);
  std::iota(indices.begin(), indices.end(), 0);
  lambdasTensor->write(indices.size(), indices.data(), lambdas.data());
  allocatedTensorArgument<F>(name + "EigenValues", lambdasTensor);
}

template <typename F>
void NaturalTransitionOrbitalsFromRhoAI::run() {
  auto RhoAI(getTensorArgument<F>("RhoAI"));
  int No(RhoAI->lens[1]), Nv(RhoAI->lens[0]);
  int oo[] = {No, No}, vv[] = {Nv, Nv}, syms[] = {NS, NS};

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "No: " << No << std::endl;
  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "Nv: " << Nv << std::endl;

  auto I = NEW(Tensor<F>, 2, oo, syms, *Sisi4s::world, "I");
  auto A = NEW(Tensor<F>, 2, vv, syms, *Sisi4s::world, "A");

  auto RhoAIConj = NEW(Tensor<F>, RhoAI);
  conjugate(*RhoAIConj);

  (*I)["ij"] = (*RhoAIConj)["ei"] * (*RhoAI)["ej"];
  (*A)["ab"] = (*RhoAIConj)["am"] * (*RhoAI)["bm"];

  buildTransformations(*I, "Occupied");
  buildTransformations(*A, "Virtual");
}
