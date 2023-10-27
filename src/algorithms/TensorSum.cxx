#include <algorithms/TensorSum.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_ALGORITHM(TensorSum) {
  Tensor<double> *A(getTensorArgument("A"));
  Tensor<double> *B(getTensorArgument("B"));
  Tensor<double> *C(getTensorArgument("Result"));
  (*C)[getTextArgument("ResultIndex").c_str()] =
      getRealArgument("AFactor") * (*A)[getTextArgument("AIndex").c_str()]
      + getRealArgument("BFactor") * (*B)[getTextArgument("BIndex").c_str()];
}
