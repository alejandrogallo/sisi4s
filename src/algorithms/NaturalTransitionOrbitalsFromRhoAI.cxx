#include <algorithms/NaturalTransitionOrbitalsFromRhoAI.hpp>
#include <util/SharedPointer.hpp>
#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <util/Log.hpp>
#include <vector>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <numeric>      // std::iota
#include <algorithm>

#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(NaturalTransitionOrbitalsFromRhoAI);

void
NaturalTransitionOrbitalsFromRhoAI::run() {
  if (getIntegerArgument("complex", 1) == 1) {
    run<cc4s::complex>();
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
      norm += std::conj(A(i,j)) * A(i,j);
    }
    LOG(1, "NaturalTransitionOrbitalsFromRhoAI") <<
      "|" << name << "(" << j << ")|**2 " << norm << std::endl;
  }
}

template <typename F>
inline void logOverlap(LapackMatrix<F> &A, const std::string name) {
  F o;
  for (int j(0); j < A.getColumns(); j++) {
    for (int i(0); i < A.getColumns(); i++) {
      o = F(0);
      for (int k(0); k < A.getRows(); k++) {
        o += std::conj(A(k,i)) * A(k,j);
      }
      LOG(1, "NaturalTransitionOrbitalsFromRhoAI") <<
        "Overlap_" << name << "(" << i << ","<<j<<") " << o << std::endl;
    }
  }
}

template <typename F>
void cleanupSpinStates(CTF::Tensor<F> &t){
  int64_t nValues;
  int64_t *globalIndices;
  F *values;
  int order(t.order);

  // read the values of the tensor in the current processor
  t.read_local(&nValues, &globalIndices, &values);
  LOG(1, "NaturalTransitionOrbitalsFromRhoAI") <<
    "local values = " << nValues << std::endl;

  for (int i=0; i<nValues; i++) {
    int g = globalIndices[i];
    // global index "carry"
    int gc = g;
    // this is a selector for the case where crossterms appear
    bool up(true), down(true);
    for (int o(0); o < order; o++) {
      int modulizer = t.lens[o];
      int ijk = gc % modulizer;
      up = up && (ijk % 2 == 0);
      down = down && ((ijk % 2 + 1) == 1);
      gc /= modulizer;
    }
    // if up and down are zero, it means that at some point the the indices
    // were part of a cross-term element, this should the be set to zero
    // because we're filtering out the terms of this kind.
    if (!(up || down)) { values[i] = F(0); }
  }

  // now we have to write the values back into the tensor
  LOG(1, "NaturalTransitionOrbitalsFromRhoAI") <<
    "Writing the cleaned up values back into the tensor" << std::endl;
  t.write(nValues, globalIndices, values);

  // clean up the mess
  delete[] values;
  delete[] globalIndices;
}

template <typename F>
void
NaturalTransitionOrbitalsFromRhoAI::buildTransformations(CTF::Tensor<F> &rho, const std::string name) {
  // We build first a lapack matrix in order to do the diagonalization
  //
  LapackMatrix<F> rhoMatrix(rho);
  int n(rho.lens[0]);
  int nn[] = {n,n};
  std::vector<int64_t> indices;

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "Diagonalizing " <<
    name << " " << n << "x" << n << " block" << std::endl;
  LapackGeneralEigenSystem<F> solver(rhoMatrix);

  // get the right eigenvectors
  LapackMatrix<complex> rightEigenVectors(solver.getRightEigenVectors());
  CTF::Tensor<F> *rightEigenVectorsTensor =
    new CTF::Tensor<F>(2, nn, rho.sym, *Cc4s::world, "r");
  // the indices will hold the vector indices
  if (rightEigenVectorsTensor->wrld->rank == 0) {
    indices.resize(n * n);
  } else {
    indices.resize(0);
  }
  std::iota(indices.begin(), indices.end(), 0);
  // write the data in the rightEigenVectors into the tensor
  rightEigenVectorsTensor->write(
    indices.size(),
    indices.data(),
    rightEigenVectors.getValues());
  allocatedTensorArgument<F>(
    name + "TransformationMatrix", rightEigenVectorsTensor);
  // log vector norms and overlaps
  logVectorNorms<F>(rightEigenVectors, name);
  logOverlap<F>(rightEigenVectors, name);

  std::vector<complex> lambdas(solver.getEigenValues());
  CTF::Tensor<F> *lambdasTensor =
    new CTF::Tensor<F>(1, nn, rho.sym, *Cc4s::world, "lambdas");
  indices.resize(n);
  std::iota(indices.begin(), indices.end(), 0);
  lambdasTensor->write(indices.size(), indices.data(), lambdas.data());
  allocatedTensorArgument<F>(name + "EigenValues", lambdasTensor);

}

template <typename F> void
NaturalTransitionOrbitalsFromRhoAI::run() {
  bool cleanup(getIntegerArgument("cleanupSpinChannels", 0) == 1);
  auto RhoAI(getTensorArgument<F>("RhoAI"));
  int No(RhoAI->lens[1]), Nv(RhoAI->lens[0]);
  int oo[] = {No, No}, vv[] = {Nv, Nv}, syms[] = {NS, NS};

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "No: " << No << std::endl;
  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "Nv: " << Nv << std::endl;

  if (cleanup) {
    LOG(0, "NaturalTransitionOrbitalsFromRhoAI") <<
      "Cleaning the Spin contamination channels" << std::endl;
    cleanupSpinStates(*RhoAI);
  }

  auto I = NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "I");
  auto A = NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "A");

  auto RhoAIConj = NEW(CTF::Tensor<F>, RhoAI);
  conjugate(*RhoAIConj);

  (*I)["ij"] = (*RhoAIConj)["ei"] * (*RhoAI)["ej"];
  (*A)["ab"] = (*RhoAIConj)["am"] * (*RhoAI)["bm"];

  buildTransformations(*I, "Occupied");
  buildTransformations(*A, "Virtual");

}
