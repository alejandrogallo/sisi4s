#include <algorithms/PseudoInverse.hpp>
#include <math/PseudoInverseHermitianSvd.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(PseudoInverse) {}


DEFSPEC(PseudoInverse,
        SPEC_IN({"A", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT({"InverseA", SPEC_VAROUT("TODO: DOC", Tensor<complex> *)}));

IMPLEMENT_ALGORITHM(PseudoInverse) {
  Tensor<complex> *A(in.get<Tensor<complex> *>("A"));
  // FIXME: use cast operators provided by CTF as soon as supported
  if (A->order != 2) throw new EXCEPTION("Matrix expected as argument A");
  CTF::Matrix<complex> *MatrixA(static_cast<CTF::Matrix<complex> *>(A));
  PseudoInverseHermitianSvd<complex> pseudoInverse(*MatrixA);
  Tensor<complex> *inverseA(new Tensor<complex>(pseudoInverse.get()));
  out.set<Tensor<complex> *>("InverseA", inverseA);
}
