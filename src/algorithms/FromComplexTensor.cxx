#include <algorithms/FromComplexTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

/**
 * \brief Testing environement
 */

DEFSPEC(FromComplexTensor,
        SPEC_IN({"A", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT({"imagA", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
                 {"RealA", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(FromComplexTensor) {
  Tensor<complex> *A(in.get<Tensor<complex> *>("A"));
  Tensor<double> *RealA(
      new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "RealA"));
  if (isArgumentGiven("imagA")) {
    Tensor<double> *ImagA(
        new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "ImagA"));
    fromComplexTensor(*A, *RealA, *ImagA);
    out.set<Tensor<double> *>("imagA", ImagA);
  } else {
    fromComplexTensor(*A, *RealA);
  }
  out.set<Tensor<double> *>("RealA", RealA);
}
