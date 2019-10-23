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

template <typename F> void
NaturalTransitionOrbitalsFromRhoAI::run() {
  auto RhoAI(getTensorArgument<F>("RhoAI"));
  int No(RhoAI->lens[1]), Nv(RhoAI->lens[0]);
  int oo[] = {No, No}, vv[] = {Nv, Nv}, syms[] = {NS, NS};

  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "No: " << No << std::endl;
  LOG(0, "NaturalTransitionOrbitalsFromRhoAI") << "Nv: " << Nv << std::endl;

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
  indices.resize(No * No);
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
  indices.resize(Nv * Nv);
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
