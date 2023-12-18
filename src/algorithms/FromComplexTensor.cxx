#include <Step.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

/**
 * \brief Testing environement
 */

DEFSPEC(FromComplexTensor,
        SPEC_IN({"data", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
        SPEC_OUT({"imag", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
                 {"real", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

DEFSTEP(FromComplexTensor) {
  Tensor<sisi4s::complex> *A(in.get<Tensor<sisi4s::complex> *>("data"));
  Tensor<double> *RealA(
      new Tensor<double>(A->order, A->lens, A->sym, *A->wrld, "RealA"));
  if (out.present("imag")) {
    Tensor<double> *ImagA(
        new Tensor<double>(A->order, A->lens, A->sym, *A->wrld, "ImagA"));
    fromComplexTensor(*A, *RealA, *ImagA);
    out.set<Tensor<double> *>("imag", ImagA);
  } else {
    fromComplexTensor(*A, *RealA);
  }
  out.set<Tensor<double> *>("real", RealA);
}
