#include <algorithms/ComplexTensorNorm.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(ComplexTensorNorm) {}


DEFSPEC(ComplexTensorNorm,
        SPEC_IN({"A", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ComplexTensorNorm) {
  Tensor<complex> *A(in.get<Tensor<complex> *>("A"));
  double norm(frobeniusNorm(*A));
  LOG(0) << "|A| = " << norm << std::endl;
  setRealArgument("Norm", norm);
}
