#include <algorithms/MeanCorrelationHoleDepth.hpp>
#include <string>
#include <vector>
#include <math/MathFunctions.hpp>
#include <algorithm>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Integrals.hpp>
#include <util/Emitter.hpp>
#include <iostream>
#include <util/Tensor.hpp>
#include <numeric>
#define LOGGER(_l) LOG(_l, "MeanCorrelationHoleDepth")

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(MeanCorrelationHoleDepth);

IMPLEMENT_EMPTY_DRYRUN(MeanCorrelationHoleDepth) {}

void MeanCorrelationHoleDepth::run() {

  checkArgumentsOrDie(
      {"PPHHDelta", "SinglesAmplitudes", "DoublesAmplitudes", "G"});

  const auto d(getTensorArgument<double>("PPHHDelta"));
  const auto T(getTensorArgument<double>("DoublesAmplitudes"));
  int No(d->lens[2]), Nv(d->lens[0]);
  std::vector<int> oo({No, No});
  std::vector<int> syms({NS, NS});
  auto gij(new Tensor<double>(2, oo.data(), syms.data(), *Sisi4s::world));

  LOGGER(0) << "Exporting G" << std::endl;
  LOGGER(0) << "No: " << No << std::endl;
  LOGGER(0) << "Nv: " << Nv << std::endl;

  EMIT() << YAML::Key << "No" << YAML::Value << No << YAML::Key << "Nv"
         << YAML::Value << Nv;

  (*gij)["ij"] = (*d)["abij"] * (*T)["abij"];

  if (isArgumentGiven("SinglesAmplitudes")) {
    const auto t(getTensorArgument<double>("SinglesAmplitudes"));
    LOGGER(0) << "Using singles amplitudes" << std::endl;
    (*gij)["ij"] += (*t)["ai"] * (*t)["bi"] * (*d)["abij"];
  }

  allocatedTensorArgument<double>("G", gij);
}
