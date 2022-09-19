#include <algorithms/GenerateRandomComplexMatrix.hpp>
#include <math/RandomTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


ALGORITHM_REGISTRAR_DEFINITION(GenerateRandomComplexMatrix);

GenerateRandomComplexMatrix::GenerateRandomComplexMatrix(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

GenerateRandomComplexMatrix::~GenerateRandomComplexMatrix() {
}

/**
 * \brief Testing environement
 */
void GenerateRandomComplexMatrix::run() {
  int m(getIntegerArgument("m", 200)), n(getIntegerArgument("n", 200));
  std::string symmetry(getTextArgument("symmetric", "none"));
  int sym(NS);
  if (symmetry == "hermitian") {
    throw new EXCEPTION("Hermitian symmetry of complex tensors not yet supported.");
  }
  Matrix<complex> *C(new Matrix<complex>(m, n, sym, *Sisi4s::world, "C"));
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  allocatedTensorArgument<complex>("Result", C);
}

