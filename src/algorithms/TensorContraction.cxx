#include <algorithms/TensorContraction.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorContraction);

TensorContraction::TensorContraction(std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

TensorContraction::~TensorContraction() {}

/**
 * \brief Testing environement
 */
void TensorContraction::run() {
  Tensor<double> *A(getTensorArgument<>("A"));
  Tensor<double> *B(getTensorArgument<>("B"));
  Tensor<double> *C(getTensorArgument<>("Result"));
  C->contract(getRealArgument("alpha", 1.0),
              *A,
              getTextArgument("AIndex").c_str(),
              *B,
              getTextArgument("BIndex").c_str(),
              getRealArgument("beta", 0.0),
              getTextArgument("ResultIndex").c_str());
}
