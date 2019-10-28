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
inline void logVectorNorms(LapackMatrix<F> &A, const char *name) {
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
inline void logOverlap(LapackMatrix<F> &A, const char *name) {
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
    if (!(up || down)) { values[g] = F(0); }
  }

  // now we have to write the values back into the tensor
  t.write(nValues, globalIndices, values);

  // clean up the mess
  delete[] values;
  delete[] globalIndices;
}

template <typename F> void
NaturalTransitionOrbitalsFromRhoAI::run() {
  auto cleanup(getIntegerArgument("cleanup-spin-channels", 0) == 1);
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

  LapackMatrix<F> IMatrix(*I);
  LapackMatrix<F> AMatrix(*A);

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") <<
    "Diagonalizing occupied " << No << "x" << No << " block" << std::endl;
  LapackGeneralEigenSystem<F> isolver(IMatrix);

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") <<
    "Diagonalizing virtual  " << Nv << "x" << Nv << " block" << std::endl;
  LapackGeneralEigenSystem<F> asolver(AMatrix);

  std::vector<int64_t> indices;

  LapackMatrix<complex> iRightEigenVectors(isolver.getRightEigenVectors());
  CTF::Tensor<F> *iRightEigenVectorsTensor = new CTF::Tensor<F>(
    2, oo, syms, *Cc4s::world, "RightEigenVectorsOccupied");
  if (iRightEigenVectorsTensor->wrld->rank == 0) {
    indices.resize(No * No);
  } else {
    indices.resize(0);
  }
  std::iota(indices.begin(), indices.end(), 0);
  iRightEigenVectorsTensor->write(
    indices.size(),
    indices.data(),
    iRightEigenVectors.getValues()
  );
  allocatedTensorArgument<F>(
    "OccupiedTransformationMatrix", iRightEigenVectorsTensor);
  logVectorNorms<F>(iRightEigenVectors, "occ");
  logOverlap<F>(iRightEigenVectors, "occ");

  LapackMatrix<complex> aRightEigenVectors(asolver.getRightEigenVectors());
  CTF::Tensor<F> *aRightEigenVectorsTensor = new CTF::Tensor<F>(
    2, vv, syms, *Cc4s::world, "RightEigenVectorsVirtual");
  if (iRightEigenVectorsTensor->wrld->rank == 0) {
    indices.resize(Nv * Nv);
  } else {
    indices.resize(0);
  }
  std::iota(indices.begin(), indices.end(), 0);
  aRightEigenVectorsTensor->write(
    indices.size(),
    indices.data(),
    aRightEigenVectors.getValues()
  );
  allocatedTensorArgument<F>(
    "VirtualTransformationMatrix", aRightEigenVectorsTensor);
  logVectorNorms<F>(aRightEigenVectors, "vir");
  logOverlap<F>(aRightEigenVectors, "vir");

  std::vector<complex> iLambdas(isolver.getEigenValues());
  CTF::Tensor<F> *iLambdasTensor = new CTF::Tensor<F>(
    1, oo, syms, *Cc4s::world, "lambdas");
  indices.resize(No);
  std::iota(indices.begin(), indices.end(), 0);
  iLambdasTensor->write(indices.size(), indices.data(), iLambdas.data());
  allocatedTensorArgument<F>("OccupiedEigenValues", iLambdasTensor);

  std::vector<complex> aLambdas(asolver.getEigenValues());
  CTF::Tensor<F> *aLambdasTensor = new CTF::Tensor<F>(
    1, vv, syms, *Cc4s::world, "lambdas");
  indices.resize(Nv);
  std::iota(indices.begin(), indices.end(), 0);
  aLambdasTensor->write(indices.size(), indices.data(), aLambdas.data());
  allocatedTensorArgument<F>("VirtualEigenValues", aLambdasTensor);



}
