#include <Step.hpp>
#include <string>
#include <vector>
#include <math/MathFunctions.hpp>
#include <algorithm>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <util/Emitter.hpp>
#include <util/SharedPointer.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <util/Tensor.hpp>
#include <numeric>
#include <map>
#include <set>
#define LOGGER(_l) LOG(_l, "FockMatrixFromCoulombIntegrals")

using namespace sisi4s;

DEFSPEC(
    FockMatrixFromCoulombIntegrals,
    SPEC_IN({"HHHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"HHMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PHHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PHMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PHPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PPMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT({"HHFockMatrix", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"HoleEigenEnergies", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"HPFockMatrix", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"ParticleEigenEnergies",
              SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"PHFockMatrix", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"PPFockMatrix", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

DEFSTEP(FockMatrixFromCoulombIntegrals) {

  const auto phph(in.get<Tensor<double> *>("PHPHCoulombIntegrals"));
  int No(phph->lens[1]), Nv(phph->lens[0]);
  std::vector<int> vv({Nv, Nv}), oo({No, No}), ov({No, Nv}), vo({Nv, No});
  std::vector<int> syms({NS, NS});

  LOGGER(0) << "No: " << No << std::endl;
  LOGGER(0) << "Nv: " << Nv << std::endl;

  // ij: HH HHHH
  auto fij(new Tensor<double>(2, oo.data(), syms.data(), *Sisi4s::world));
  const auto hh(in.get<Tensor<double> *>("HHMatrix"));
  const auto hhhh(in.get<Tensor<double> *>("HHHHCoulombIntegrals"));
  (*fij)["ij"] = (*hh)["ij"];
  (*fij)["ij"] += (+2.0) * (*hhhh)["ikjk"];
  (*fij)["ij"] += (-1.0) * (*hhhh)["ikkj"];
  out.set<Tensor<double> *>("HHFockMatrix", fij);

  // ab: PP PHPH PHHP
  auto fab(new Tensor<double>(2, vv.data(), syms.data(), *Sisi4s::world));
  const auto pp(in.get<Tensor<double> *>("PPMatrix"));
  const auto phhp(in.get<Tensor<double> *>("PHHPCoulombIntegrals"));
  (*fab)["ab"] = (*pp)["ab"];
  (*fab)["ab"] += (+2.0) * (*phph)["akbk"];
  (*fab)["ab"] += (-1.0) * (*phhp)["akkb"];
  out.set<Tensor<double> *>("PPFockMatrix", fab);

  // ai: PH PHHH
  auto fai(new Tensor<double>(2, vo.data(), syms.data(), *Sisi4s::world));
  const auto ph(in.get<Tensor<double> *>("PHMatrix"));
  const auto phhh(in.get<Tensor<double> *>("PHHHCoulombIntegrals"));
  (*fai)["ai"] = (*ph)["ai"];
  (*fai)["ai"] += (+2.0) * (*phhh)["akik"];
  (*fai)["ai"] += (-1.0) * (*phhh)["akki"];
  out.set<Tensor<double> *>("PHFockMatrix", fai);

  // ia: HP HHPH HHHP
  auto fia(new Tensor<double>(2, ov.data(), syms.data(), *Sisi4s::world));
  (*fia)["ia"] = (*fai)["ai"];
  out.set<Tensor<double> *>("HPFockMatrix", fia);

  auto epsi(new Tensor<double>(1, oo.data(), syms.data(), *Sisi4s::world));
  (*epsi)["i"] = (*fij)["ii"];
  out.set<Tensor<double> *>("HoleEigenEnergies", epsi);

  auto epsa(new Tensor<double>(1, vv.data(), syms.data(), *Sisi4s::world));
  (*epsa)["a"] = (*fab)["aa"];
  out.set<Tensor<double> *>("ParticleEigenEnergies", epsa);

  CTF::Scalar<double> energy;
  energy[""] = (+2.0) * (*hh)["ii"];
  energy[""] += (+2.0) * (*hhhh)["ikik"];
  energy[""] += (-1.0) * (*hhhh)["ikki"];
  const double dEnergy(energy.get_val());

  EMIT() << YAML::Key << "energy" << YAML::Value << dEnergy << YAML::Key << "No"
         << YAML::Value << No << YAML::Key << "Nv" << YAML::Value << Nv;

  LOG(0, "FockMatrixFromCoulombIntegrals")
      << "energy= " << dEnergy << std::endl;
}
