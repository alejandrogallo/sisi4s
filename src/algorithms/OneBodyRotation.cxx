#include <Step.hpp>

using namespace sisi4s;
#define LOGGER(_l) LOG(_l, "OneBodyRotation")

// TODO: Update spec for compelx and double
DEFSPEC(
    OneBodyRotation,
    SPEC_IN({"No",
             SPEC_POSITIVE("The number of occupied orbitals", int)->require()},
            {"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},
            {"OrbitalCoefficients", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT({"hh",
              SPEC_VAROUT("Get the hole-hole part of the operator",
                          Tensor<double> *)},
             {"pp",
              SPEC_VAROUT("Get the particle-particle part of the operator",
                          Tensor<double> *)},
             {"hp",
              SPEC_VAROUT("Get the hole-particle part of the operator",
                          Tensor<double> *)},
             {"ph",
              SPEC_VAROUT("Get the particle-hole part of the operator",
                          Tensor<double> *)},
             {"Out",
              SPEC_VAROUT("Do not slice, this is the whole operator",
                          Tensor<double> *)}));

template <typename F>
static void run(Arguments &in, Arguments &out) {
  auto C(in.get<Tensor<F> *>("OrbitalCoefficients"));
  auto I(in.get<Tensor<F> *>("Data"));
  auto O(new Tensor<F>(*I));
  const int No(in.get<int>("No")), Np(C->lens[0]);

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
    if (!out.present(i.name)) continue;
    auto tensor(new Tensor<F>(O->slice(i.begin, i.end)));
    LOGGER(0) << i.name << ": "
              << "{" << i.begin[0] << "," << i.begin[1] << "}"
              << " --> "
              << "{" << i.end[0] << "," << i.end[1] << "}" << std::endl;
    out.set<Tensor<F> *>(i.name, tensor);
  }
}

DEFSTEP(OneBodyRotation) {
  if (in.is_of_type<Tensor<double> *>("OrbitalCoefficients")) {
    ::run<double>(in, out);
  } else {
    ::run<sisi4s::complex>(in, out);
  }
}
