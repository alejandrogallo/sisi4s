#include <algorithms/ComplexTensorNorm.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorNorm);

ComplexTensorNorm::ComplexTensorNorm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ComplexTensorNorm::~ComplexTensorNorm() {
}

/**
 * \brief Testing environement
 */
void ComplexTensorNorm::run() {
  Tensor<complex> *A(getTensorArgument<complex>("A"));
  double norm(frobeniusNorm(*A));
  LOG(0) << "|A| = " << norm << std::endl;
  setRealArgument("Norm", norm);
}
