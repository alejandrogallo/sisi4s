#include <algorithms/ComplexTensorContraction.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorContraction);

/**
 * \brief Testing environement
 */
void ComplexTensorContraction::run() {

  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<complex> *B(getTensorArgument<complex>("B"));
  Tensor<complex> *C(getTensorArgument<complex>("Result"));
  (*C)[getTextArgument("ResultIndex").c_str()] =
      (*A)[getTextArgument("AIndex").c_str()]
      * (*B)[getTextArgument("BIndex").c_str()];
}
