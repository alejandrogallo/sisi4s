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

IMPLEMENT_EMPTY_DRYRUN(MeanCorrelationHoleDepth) {}

DEFSPEC(
    MeanCorrelationHoleDepth,
    SPEC_IN({"DoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PPHHDelta", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"SinglesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT({"G", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(MeanCorrelationHoleDepth) {

  const auto d(in.get<Tensor<double> *>("PPHHDelta"));
  const auto T(in.get<Tensor<double> *>("DoublesAmplitudes"));
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

  if (in.present("SinglesAmplitudes")) {
    const auto t(in.get<Tensor<double> *>("SinglesAmplitudes"));
    LOGGER(0) << "Using singles amplitudes" << std::endl;
    (*gij)["ij"] += (*t)["ai"] * (*t)["bi"] * (*d)["abij"];
  }

  out.set<Tensor<double> *>("G", gij);
}
