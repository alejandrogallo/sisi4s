#include <algorithms/FromComplexTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(FromComplexTensor);

/**
 * \brief Testing environement
 */
void FromComplexTensor::run() {
  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<double> *RealA(
      new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "RealA"));
  if (isArgumentGiven("imagA")) {
    Tensor<double> *ImagA(
        new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "ImagA"));
    fromComplexTensor(*A, *RealA, *ImagA);
    allocatedTensorArgument("imagA", ImagA);
  } else {
    fromComplexTensor(*A, *RealA);
  }
  allocatedTensorArgument("RealA", RealA);
}
