#include <Step.hpp>
#include <math/PseudoInverseHermitianSvd.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(PseudoInverse,
        SPEC_IN({"A", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
        SPEC_OUT({"InverseA",
                  SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)}));

DEFSTEP(PseudoInverse) {
  Tensor<sisi4s::complex> *A(in.get<Tensor<sisi4s::complex> *>("A"));
  // FIXME: use cast operators provided by CTF as soon as supported
  if (A->order != 2) throw new EXCEPTION("Matrix expected as argument A");
  CTF::Matrix<complex> *MatrixA(static_cast<CTF::Matrix<complex> *>(A));
  PseudoInverseHermitianSvd<complex> pseudoInverse(*MatrixA);
  Tensor<sisi4s::complex> *inverseA(
      new Tensor<sisi4s::complex>(pseudoInverse.get()));
  out.set<Tensor<sisi4s::complex> *>("InverseA", inverseA);
}
