#include <algorithms/ComplexTensorSum.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorSum);

/**
 * \brief Testing environement
 */
void ComplexTensorSum::run() {

  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<complex> *B(getTensorArgument<complex>("B"));
  Tensor<complex> *C(getTensorArgument<complex>("Result"));
  (*C)[getTextArgument("ResultIndex").c_str()] =
      getRealArgument("AFactor") * (*A)[getTextArgument("AIndex").c_str()]
      + getRealArgument("BFactor") * (*B)[getTextArgument("BIndex").c_str()];
}
