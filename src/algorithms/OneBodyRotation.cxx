#include <vector>
#include <algorithms/OneBodyRotation.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(OneBodyRotation) {}

#define LOGGER(_l) LOG(_l, "OneBodyRotation")

DEFSPEC(OneBodyRotation,
        SPEC_IN({"No", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
                {"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
                {"OrbitalCoefficients",
                 SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({i.name, SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(OneBodyRotation) {

  auto C(in.get<Tensor<double> *>("OrbitalCoefficients"));
  auto I(in.get<Tensor<double> *>("Data"));
  auto O(new Tensor<double>(*I));
  const int No(in.get<int64_t>("No", 0)), Np(C->lens[0]);

  LOGGER(0) << "No: " << No << std::endl;
  LOGGER(0) << "Np: " << Np << std::endl;

  LOGGER(0) << "Rotating" << std::endl;
  (*O)["pq"] = (*C)["Pp"] * (*C)["Qq"] * (*I)["PQ"];

  struct _info {
    const char *name;
    const int begin[2], end[2];
  };
  const std::vector<_info> infos = {{"hh", {0, 0}, {No, No}},
                                    {"pp", {No, No}, {Np, Np}},
                                    {"hp", {0, No}, {No, Np}},
                                    {"ph", {No, 0}, {Np, No}},
                                    {"Out", {0, 0}, {Np, Np}}};

  for (const auto &i : infos) {
    if (!isArgumentGiven(i.name)) continue;
    auto tensor(new Tensor<double>(O->slice(i.begin, i.end)));
    LOGGER(0) << i.name << ": "
              << "{" << i.begin[0] << "," << i.begin[1] << "}"
              << " --> "
              << "{" << i.end[0] << "," << i.end[1] << "}" << std::endl;
    out.set<Tensor<double> *>(i.name, tensor);
  }
}
