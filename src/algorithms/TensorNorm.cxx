#include <algorithms/TensorNorm.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorNorm) {}

DEFSPEC(TensorNorm,
        SPEC_IN({"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"Norm", SPEC_VAROUT("TODO: DOC", double)}));

IMPLEMENT_ALGORITHM(TensorNorm) {
  Tensor<double> *A(in.get<Tensor<double> *>("Data"));
  double norm(frobeniusNorm(*A));
  LOG(0, "TensorNorm") << "norm = " << norm << std::endl;
  if (out.present("Norm")) out.set<double>("Norm", norm);
}
