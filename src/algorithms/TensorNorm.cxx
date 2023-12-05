#include <algorithms/TensorNorm.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorNorm) {}

IMPLEMENT_ALGORITHM(TensorNorm) {
  Tensor<double> *A(getTensorArgument("Data"));
  double norm(frobeniusNorm(*A));
  LOG(0, "TensorNorm") << "norm = " << norm << std::endl;
  if (isArgumentGiven("Norm")) setRealArgument("Norm", norm);
}
