#include <algorithms/MeanCorrelationHoleDepth.hpp>
#include <string>
#include <vector>
#include <math/MathFunctions.hpp>
#include <algorithm>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Integrals.hpp>
#include <util/Emitter.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#define LOGGER(_l) LOG(_l, "MeanCorrelationHoleDepth")

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(MeanCorrelationHoleDepth);

void MeanCorrelationHoleDepth::run() {

  checkArgumentsOrDie({"PPHHDelta", "DoublesAmplitudes", "G"});

  const auto d(getTensorArgument<double>("PPHHDelta"));
  const auto t(getTensorArgument<double>("DoublesAmplitudes"));
  int No(d->lens[2]), Nv(d->lens[0]);
  std::vector<int> oo({No, No});
  std::vector<int> syms({NS, NS});
  auto gij(new CTF::Tensor<double>(2, oo.data(), syms.data(), *Cc4s::world));

  LOGGER(0) << "Exporting G" << std::endl;
  LOGGER(0) << "No: " << No << std::endl;
  LOGGER(0) << "Nv: " << Nv << std::endl;

  EMIT() << YAML::Key << "No" << YAML::Value << No
         << YAML::Key << "Nv" << YAML::Value << Nv
         ;

  (*gij)["ij"] = (*d)["abij"] * (*t)["abij"];

  allocatedTensorArgument<double>("G", gij);

}
