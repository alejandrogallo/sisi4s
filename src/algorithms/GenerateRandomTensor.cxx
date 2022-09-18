#include <algorithms/GenerateRandomTensor.hpp>
#include <math/RandomTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


ALGORITHM_REGISTRAR_DEFINITION(GenerateRandomTensor);

GenerateRandomTensor::GenerateRandomTensor(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

GenerateRandomTensor::~GenerateRandomTensor() {
}

/**
 * \brief Testing environement
 */
void GenerateRandomTensor::run() {

  int nv(getIntegerArgument("Nv")), no(getIntegerArgument("No")); 
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<double> *C(new Tensor<double>(4, lens, syms, *Sisi4s::world, "C"));
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  allocatedTensorArgument("Result", C);
  
}
