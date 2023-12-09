#include <system_error>
#include <algorithms/CoulombIntegralsFromVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

using namespace sisi4s;

DEFSPEC(
    CoulombIntegralsFromVertex,
    SPEC_IN({"antisymmetrize", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
            {"complex", SPEC_VALUE_DEF("TODO: DOC", bool, 0)},
            {"forceReal", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
            {"CoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
            {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"ParticleEigenEnergies",
             SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT(
        {"HHHHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HHHPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HHPHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HHPPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HPHHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HPHPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HPPHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"HPPPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PHHHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PHHPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PHPHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PHPPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PPHHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PPHPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PPPHCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PPPPCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"PQRSCoulombIntegrals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(CoulombIntegralsFromVertex) {
  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr(in.get<Tensor<complex> *>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));
  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));

  LOG(0, "Integrals") << "Reading Coulomb integrals form vertex "
                      << GammaGqr->get_name() << std::endl;

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(GammaGqr->lens[1]);
  EMIT() << YAML::Key << "NG" << YAML::Value << NG << YAML::Key << "No"
         << YAML::Value << No << YAML::Key << "Nv" << YAML::Value << Nv
         << YAML::Key << "Np" << YAML::Value << Np;

  std::vector<std::string> list;
  if (out.present("PHPHCoulombIntegrals")) {
    list.push_back("PHPHCoulombIntegrals");
  }
  if (out.present("PPHHCoulombIntegrals")) {
    list.push_back("PPHHCoulombIntegrals");
  }
  if (out.present("HHHHCoulombIntegrals")) {
    list.push_back("HHHHCoulombIntegrals");
  }
  if (out.present("HHHPCoulombIntegrals")) {
    list.push_back("HHHPCoulombIntegrals");
  }
  if (out.present("PPPPCoulombIntegrals")) {
    list.push_back("PPPPCoulombIntegrals");
  }
  if (out.present("PPPHCoulombIntegrals")) {
    list.push_back("PPPHCoulombIntegrals");
  }
  if (out.present("PHHHCoulombIntegrals")) {
    list.push_back("PHHHCoulombIntegrals");
  }
  if (out.present("HHPPCoulombIntegrals")) {
    list.push_back("HHPPCoulombIntegrals");
  }
  if (out.present("PHHPCoulombIntegrals")) {
    list.push_back("PHHPCoulombIntegrals");
  }
  if (out.present("HPHHCoulombIntegrals")) {
    list.push_back("HPHHCoulombIntegrals");
  }
  if (out.present("HPHPCoulombIntegrals")) {
    list.push_back("HPHPCoulombIntegrals");
  }
  if (out.present("HPPPCoulombIntegrals")) {
    list.push_back("HPPPCoulombIntegrals");
  }
  if (out.present("PPHPCoulombIntegrals")) {
    list.push_back("PPHPCoulombIntegrals");
  }
  if (out.present("HPPHCoulombIntegrals")) {
    list.push_back("HPPHCoulombIntegrals");
  }
  if (out.present("HHPHCoulombIntegrals")) {
    list.push_back("HHPHCoulombIntegrals");
  }
  if (out.present("PHPPCoulombIntegrals")) {
    list.push_back("PHPPCoulombIntegrals");
  }
  EMIT() << YAML::Key << "integrals" << YAML::Value;
  EMIT() << YAML::Flow << YAML::BeginSeq;
  for (size_t i(0); i < list.size(); i++) { EMIT() << list[i]; }
  EMIT() << YAML::EndSeq;

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  // TODO: calculate vvvv, vvvo  for COMPLEX
  syms = std::array<int, 4>{{NS, NS, NS, NS}};
  vvvv = std::array<int, 4>{{Nv, Nv, Nv, Nv}};
  vovo = std::array<int, 4>{{Nv, No, Nv, No}};
  vvoo = std::array<int, 4>{{Nv, Nv, No, No}};
  oooo = std::array<int, 4>{{No, No, No, No}};
  ooov = std::array<int, 4>{{No, No, No, Nv}};
  vvvo = std::array<int, 4>{{Nv, Nv, Nv, No}};
  vooo = std::array<int, 4>{{Nv, No, No, No}};
  voov = std::array<int, 4>{{Nv, No, No, Nv}};
  oovv = std::array<int, 4>{{No, No, Nv, Nv}};
  ovoo = std::array<int, 4>{{No, Nv, No, No}};
  ovov = std::array<int, 4>{{No, Nv, No, Nv}};
  ovvv = std::array<int, 4>{{No, Nv, Nv, Nv}};
  vvov = std::array<int, 4>{{Nv, Nv, No, Nv}};
  ovvo = std::array<int, 4>{{No, Nv, Nv, No}};
  oovo = std::array<int, 4>{{No, No, Nv, No}};
  vovv = std::array<int, 4>{{Nv, No, Nv, Nv}};

  bool realIntegrals = !in.get<bool>("complex");
  LOG(0, "CoulombIntegrals") << "Using " << (realIntegrals ? "real" : "complex")
                             << " Coulomb integrals" << std::endl;

  int aStart(Np - Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int GijStart[] = {0, iStart, iStart};
  int GijEnd[] = {NG, iEnd, iEnd};
  int GiaStart[] = {0, iStart, aStart};
  int GiaEnd[] = {NG, iEnd, aEnd};
  int GaiStart[] = {0, aStart, iStart};
  int GaiEnd[] = {NG, aEnd, iEnd};
  int GabStart[] = {0, aStart, aStart};
  int GabEnd[] = {NG, aEnd, aEnd};
  GammaGij = new Tensor<complex>(GammaGqr->slice(GijStart, GijEnd));
  GammaGia = realIntegrals
               ? nullptr
               : new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd));
  GammaGai = new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd));
  GammaGab = new Tensor<complex>(GammaGqr->slice(GabStart, GabEnd));

  if (realIntegrals) {
    calculateRealIntegrals();
  } else {
    calculateComplexIntegrals();
  }

  delete GammaGij;
  delete GammaGai;
  delete GammaGab;
  if (GammaGia) { delete GammaGia; }
}

void CoulombIntegralsFromVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(in.get<DryTensor<complex> *>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(in.get<DryTensor<double> *>("HoleEigenEnergies"));
  DryTensor<> *epsa(in.get<DryTensor<double> *>("ParticleEigenEnergies"));

  LOG(0, "Integrals") << "Reading Coulomb integrals form vertex "
                      << GammaGqr->get_name() << std::endl;

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  syms = std::array<int, 4>{{NS, NS, NS, NS}};
  vvvv = std::array<int, 4>{{Nv, Nv, Nv, Nv}};
  vovo = std::array<int, 4>{{Nv, No, Nv, No}};
  vvoo = std::array<int, 4>{{Nv, Nv, No, No}};
  voov = std::array<int, 4>{{Nv, No, No, Nv}};
  oovv = std::array<int, 4>{{No, No, Nv, Nv}};
  oooo = std::array<int, 4>{{No, No, No, No}};
  ooov = std::array<int, 4>{{No, No, No, Nv}};
  vooo = std::array<int, 4>{{Nv, No, No, No}};
  vvvo = std::array<int, 4>{{Nv, Nv, Nv, No}};

  bool realIntegrals = !in.get<bool>("complex");
  LOG(0, "CoulombIntegrals") << "Using " << (realIntegrals ? "real" : "complex")
                             << " Coulomb integrals" << std::endl;

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiLens[] = {NG, Nv, No};
  int GiaLens[] = {NG, No, Nv};
  int GabLens[] = {NG, Nv, Nv};
  int GijLens[] = {NG, No, No};

  DryTensor<complex> GammaGia(3, GiaLens, syms.data());
  DryTensor<complex> GammaGai(3, GaiLens, syms.data());
  DryTensor<complex> GammaGab(3, GabLens, syms.data());
  DryTensor<complex> GammaGij(3, GijLens, syms.data());

  if (realIntegrals) {
    dryCalculateRealIntegrals();
  } else {
    dryCalculateComplexIntegrals();
  }
}

void CoulombIntegralsFromVertex::calculateRealIntegrals() {

  const bool calculate_vpqrs = out.present("PQRSCoulombIntegrals");

  int antisymmetrize(in.get<int64_t>("antisymmetrize"));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals")
        << "Calculating antisymmetrized integrals" << std::endl;
  }
  Tensor<double> *Vaibj(out.present("PHPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vovo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaibj")
                            : nullptr);
  Tensor<double> *Vabij(out.present("PPHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvoo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabij")
                            : nullptr);
  Tensor<double> *Vijkl(out.present("HHHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oooo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijkl")
                            : nullptr);
  Tensor<double> *Vijka(out.present("HHHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ooov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijka")
                            : nullptr);
  Tensor<double> *Vabcd(out.present("PPPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvvv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabcd")
                            : nullptr);
  Tensor<double> *Vabci(out.present("PPPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvvo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabci")
                            : nullptr);

  // Initialization of tensors created from already existing ones
  Tensor<double> *Vaijk(out.present("PHHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vooo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaijk")
                            : nullptr);
  Tensor<double> *Vijab(out.present("HHPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oovv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijab")
                            : nullptr);
  Tensor<double> *Vaijb(out.present("PHHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 voov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaijb")
                            : nullptr);
  Tensor<double> *Viajk(out.present("HPHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovoo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viajk")
                            : nullptr);
  Tensor<double> *Viajb(out.present("HPHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viajb")
                            : nullptr);
  Tensor<double> *Viabc(out.present("HPPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovvv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viabc")
                            : nullptr);
  Tensor<double> *Vabic(out.present("PPHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabic")
                            : nullptr);
  Tensor<double> *Viabj(out.present("HPPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovvo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viabj")
                            : nullptr);
  Tensor<double> *Vijak(out.present("HHPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oovo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijak")
                            : nullptr);
  Tensor<double> *Vaibc(out.present("PHPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vovv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijak")
                            : nullptr);

  // Split GammaGab,GammaGai,GammaGia,GammaGij into real and imaginary parts
  Tensor<double> realGammaGai(3,
                              GammaGai->lens,
                              GammaGai->sym,
                              *GammaGai->wrld,
                              "RealGammaGai");
  Tensor<double> imagGammaGai(3,
                              GammaGai->lens,
                              GammaGai->sym,
                              *GammaGai->wrld,
                              "ImagGammaGai");
  fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

  Tensor<double> realGammaGab(3,
                              GammaGab->lens,
                              GammaGab->sym,
                              *GammaGab->wrld,
                              "RealGammaGab");
  Tensor<double> imagGammaGab(3,
                              GammaGab->lens,
                              GammaGab->sym,
                              *GammaGab->wrld,
                              "ImagGammaGab");
  fromComplexTensor(*GammaGab, realGammaGab, imagGammaGab);

  Tensor<double> realGammaGij(3,
                              GammaGij->lens,
                              GammaGij->sym,
                              *GammaGij->wrld,
                              "RealGammaGij");
  Tensor<double> imagGammaGij(3,
                              GammaGij->lens,
                              GammaGij->sym,
                              *GammaGij->wrld,
                              "ImagGammaGij");
  fromComplexTensor(*GammaGij, realGammaGij, imagGammaGij);

  if (calculate_vpqrs) {
    Tensor<complex> *GammaGqr(in.get<Tensor<complex> *>("CoulombVertex"));
    int64_t No = (GammaGij->lens[1]), Nv = (GammaGab->lens[1]);
    LOG(1, "CoulombIntegrals") << "Evaluating full Vpqrs" << std::endl;
    const auto pppp = std::vector<int64_t>(4, GammaGqr->lens[1]);
    Tensor<double> *t = new Tensor<double>(pppp.size(),
                                           pppp.data(),
                                           syms.data(),
                                           *Sisi4s::world,
                                           "Vpqrs");
    Tensor<double> realGammaGqr(3,
                                GammaGqr->lens,
                                GammaGqr->sym,
                                *GammaGqr->wrld,
                                "RealGammaGqr");
    Tensor<double> imagGammaGqr(3,
                                GammaGqr->lens,
                                GammaGqr->sym,
                                *GammaGqr->wrld,
                                "ImagGammaGqr");
    fromComplexTensor(*GammaGqr, realGammaGqr, imagGammaGqr);
    (*t)["pqrs"] = realGammaGqr["Gpr"] * realGammaGqr["Gqs"];
    (*t)["pqrs"] += imagGammaGqr["Gpr"] * imagGammaGqr["Gqs"];
    if (antisymmetrize) { (*t)["pqrs"] -= (*t)["pqsr"]; }
    out.set<Tensor<double> *>("PQRSCoulombIntegrals", t);
  }

  // Compute the integrals Vabij Vaibj Vaijb Vijkl Vabcd
  if (Vaibj) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"] = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
    out.set<Tensor<double> *>("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabij->get_name() << std::endl;
    (*Vabij)["abij"] = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
    out.set<Tensor<double> *>("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijkl) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"] = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
    out.set<Tensor<double> *>("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"] = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
    out.set<Tensor<double> *>("HHHPCoulombIntegrals", Vijka);
  }
  if (Vabcd) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"] = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
    out.set<Tensor<double> *>("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vabci) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabci->get_name() << std::endl;
    (*Vabci)["abci"] = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
    out.set<Tensor<double> *>("PPPHCoulombIntegrals", Vabci);
  }

  // Create the rest of integrals from the already given ones
  // --------------------------------------------------------

  if (Viabj) {
    // ovvo = hl * vvoo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viabj->get_name()
                               << " using " << Vabij->get_name() << std::endl;

    (*Viabj)["iabj"] = (*Vabij)["baij"];
    if (antisymmetrize) {
      // ovov = v * vovo
      (*Viabj)["iabj"] -= (*Vaibj)["aibj"];
    }
    out.set<Tensor<double> *>("HPPHCoulombIntegrals", Viabj);
  }
  if (Viajb) {
    // ovov = v * vovo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viajb->get_name()
                               << " using " << Vaibj->get_name() << std::endl;

    (*Viajb)["iajb"] = (*Vaibj)["aibj"];
    if (antisymmetrize) {
      // ovvo = hl * vvoo
      (*Viajb)["iajb"] -= (*Vabij)["baij"];
    }
    out.set<Tensor<double> *>("HPHPCoulombIntegrals", Viajb);
  }
  if (Viabc) {
    // ovvv = h%v * vvvo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viabc->get_name()
                               << " using " << Vabci->get_name() << std::endl;

    (*Viabc)["iabc"] = (*Vabci)["cbai"];
    if (antisymmetrize) {
      // ovvv = h%v * vvvo
      (*Viabc)["iabc"] -= (*Vabci)["bcai"];
    }
    out.set<Tensor<double> *>("HPPPCoulombIntegrals", Viabc);
  }
  if (Vijak) {
    // oovo = v * ooov
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vijak->get_name()
                               << " using " << Vijka->get_name() << std::endl;

    (*Vijak)["ijak"] = (*Vijka)["jika"];
    if (antisymmetrize) {
      // ooov = e * ooov
      (*Vijak)["ijak"] -= (*Vijka)["ijka"];
    }
    out.set<Tensor<double> *>("HHPHCoulombIntegrals", Vijak);
  }
  if (Vijab) {
    // oovv = h * vvoo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vijab->get_name()
                               << " using " << Vabij->get_name() << std::endl;

    (*Vijab)["ijab"] = (*Vabij)["abij"];
    if (antisymmetrize) {
      // oovv = h * vvoo
      (*Vijab)["ijab"] -= (*Vabij)["baij"];
    }
    out.set<Tensor<double> *>("HHPPCoulombIntegrals", Vijab);
  }
  if (Vabic) {
    // vvov = v * vvvo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vabic->get_name()
                               << " using " << Vabci->get_name() << std::endl;

    (*Vabic)["abic"] = (*Vabci)["baci"];
    if (antisymmetrize) {
      // vvvo = e * vvvo
      (*Vabic)["abic"] -= (*Vabci)["abci"];
    }
    out.set<Tensor<double> *>("PPHPCoulombIntegrals", Vabic);
  }
  if (Vaijb) {
    // voov = hr * vvoo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vaijb->get_name()
                               << " using " << Vabij->get_name() << std::endl;

    (*Vaijb)["aijb"] = (*Vabij)["abji"];
    if (antisymmetrize) {
      // vovo = e * vovo
      (*Vaijb)["aijb"] -= (*Vaibj)["aibj"];
    }
    out.set<Tensor<double> *>("PHHPCoulombIntegrals", Vaijb);
  }
  if (Vaibc) {
    // vovv = h * vvvo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vaibc->get_name()
                               << " using " << Vabci->get_name() << std::endl;

    (*Vaibc)["aibc"] = (*Vabci)["bcai"];
    if (antisymmetrize) {
      // vovv = h * vvvo
      (*Vaibc)["aibc"] -= (*Vabci)["cbai"];
    }
    out.set<Tensor<double> *>("PHPPCoulombIntegrals", Vaibc);
  }
  if (Vaijk) {
    // vooo = h%v * ooov
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vaijk->get_name()
                               << " using " << Vijka->get_name() << std::endl;

    (*Vaijk)["aijk"] = (*Vijka)["kjia"];
    if (antisymmetrize) {
      // vooo = h%v * ooov
      (*Vaijk)["aijk"] -= (*Vijka)["jkia"];
    }
    out.set<Tensor<double> *>("PHHHCoulombIntegrals", Vaijk);
  }
  if (Viajk) {
    // ovoo = h * ooov
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viajk->get_name()
                               << " using " << Vijka->get_name() << std::endl;

    (*Viajk)["iajk"] = (*Vijka)["jkia"];
    if (antisymmetrize) {
      // ovoo = h * ooov
      (*Viajk)["iajk"] -= (*Vijka)["kjia"];
    }
    out.set<Tensor<double> *>("HPHHCoulombIntegrals", Viajk);
  }

  if (antisymmetrize) {
    // Antisymmetrize integrals calculated directly from Gamma
    // IMPORTANT: This must be written after the creation of the integrals
    //            depending on them.
    if (Vijkl) (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
    if (Vabcd) (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
    if (Vijka) (*Vijka)["ijka"] -= (*Vijka)["jika"];
    if (Vabci) (*Vabci)["abci"] -= (*Vabci)["baci"];
    // Vaibj depends on Vabij, do not change the order
    if (Vaibj) (*Vaibj)["aibj"] -= (*Vabij)["baij"];
    if (Vabij) (*Vabij)["abij"] -= (*Vabij)["abji"];
  }
}

void CoulombIntegralsFromVertex::dryCalculateRealIntegrals() {
  DryTensor<> *Vaibj(out.present("PHPHCoulombIntegrals")
                         ? new DryTensor<>(4, vovo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabij(out.present("PPHHCoulombIntegrals")
                         ? new DryTensor<>(4, vvoo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijab(out.present("HHPPCoulombIntegrals")
                         ? new DryTensor<>(4, oovv.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vaijb(out.present("PHHPCoulombIntegrals")
                         ? new DryTensor<>(4, voov.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijkl(out.present("HHHHCoulombIntegrals")
                         ? new DryTensor<>(4, oooo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijka(out.present("HHHPCoulombIntegrals")
                         ? new DryTensor<>(4, ooov.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vaijk(out.present("PHHHCoulombIntegrals")
                         ? new DryTensor<>(4, vooo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabcd(out.present("PPPPCoulombIntegrals")
                         ? new DryTensor<>(4, vvvv.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabci(out.present("PPPHCoulombIntegrals")
                         ? new DryTensor<>(4, vvvo.data(), syms.data())
                         : nullptr);

  if (Vaibj) { out.set<DryTensor<double> *>("PHPHCoulombIntegrals", Vaibj); }
  if (Vabij) { out.set<DryTensor<double> *>("PPHHCoulombIntegrals", Vabij); }
  if (Vijab) { out.set<DryTensor<double> *>("HHPPCoulombIntegrals", Vijab); }
  if (Vaijb) { out.set<DryTensor<double> *>("PHHPCoulombIntegrals", Vaijb); }
  if (Vijkl) { out.set<DryTensor<double> *>("HHHHCoulombIntegrals", Vijkl); }
  if (Vijka) { out.set<DryTensor<double> *>("HHHPCoulombIntegrals", Vijka); }
  if (Vaijk) { out.set<DryTensor<double> *>("PHHHCoulombIntegrals", Vaijk); }
  if (Vabcd) { out.set<DryTensor<double> *>("PPPPCoulombIntegrals", Vabcd); }
  if (Vabci) { out.set<DryTensor<double> *>("PPPHCoulombIntegrals", Vabci); }
}

void CoulombIntegralsFromVertex::calculateComplexIntegrals() {
  int antisymmetrize(in.get<int64_t>("antisymmetrize"));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals")
        << "Calculating antisymmetrized integrals" << std::endl;
  }

  Tensor<complex> *Vabij(out.present("PPHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvoo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabij")
                             : nullptr);

  Tensor<complex> *Vijab(
      // TODO: HHPP is always conj(Permute(PPHH))
      out.present("HHPPCoulombIntegrals") ? new Tensor<complex>(4,
                                                                oovv.data(),
                                                                syms.data(),
                                                                *Sisi4s::world,
                                                                "Vijab")
                                          : nullptr);

  Tensor<complex> *Vaijb(out.present("PHHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   voov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaijb")
                             : nullptr);

  Tensor<complex> *Vaibj(out.present("PHPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vovo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaibj")
                             : nullptr);

  Tensor<complex> *Vijkl(out.present("HHHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   oooo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijkl")
                             : nullptr);

  Tensor<complex> *Vijka(out.present("HHHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ooov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijka")
                             : nullptr);

  Tensor<complex> *Vaijk(
      // TODO: PHHH is always conj(Permute(HHHP))
      out.present("PHHHCoulombIntegrals") ? new Tensor<complex>(4,
                                                                vooo.data(),
                                                                syms.data(),
                                                                *Sisi4s::world,
                                                                "Vaijk")
                                          : nullptr);

  Tensor<complex> *Vabcd(out.present("PPPPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvvv.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabcd")
                             : nullptr);

  Tensor<complex> *Vaibc(out.present("PHPPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vovv.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaibc")
                             : nullptr);

  // Initialization of tensors created from already existing ones
  // In principle the integrals above do not constitute the minimal set of
  // integrals from which one can write the rest (which would be 7)

  Tensor<complex> *Vabic(out.present("PPHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabic")
                             : nullptr);
  Tensor<complex> *Vabci(out.present("PPPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvvo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabci")
                             : nullptr);
  Tensor<complex> *Vijak(out.present("HHPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   oovo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijak")
                             : nullptr);
  Tensor<complex> *Viajk(out.present("HPHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovoo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viajk")
                             : nullptr);
  Tensor<complex> *Viajb(out.present("HPHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viajb")
                             : nullptr);
  Tensor<complex> *Viabj(out.present("HPPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovvo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viabj")
                             : nullptr);
  Tensor<complex> *Viabc(out.present("HPPPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovvv.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viabc")
                             : nullptr);

  CTF::Univar_Function<complex> fConj(conj<complex>);

  Tensor<complex> conjTransposeGammaGai(false, *GammaGai);
  conjTransposeGammaGai.sum(1.0, *GammaGia, "Gia", 0.0, "Gai", fConj);

  Tensor<complex> conjTransposeGammaGia(false, *GammaGia);
  conjTransposeGammaGia.sum(1.0, *GammaGai, "Gai", 0.0, "Gia", fConj);

  Tensor<complex> conjTransposeGammaGij(false, *GammaGij);
  conjTransposeGammaGij.sum(1.0, *GammaGij, "Gji", 0.0, "Gij", fConj);

  Tensor<complex> conjTransposeGammaGab(false, *GammaGab);
  conjTransposeGammaGab.sum(1.0, *GammaGab, "Gba", 0.0, "Gab", fConj);

  if (Vabij) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabij->get_name() << std::endl;
    (*Vabij)["abij"] = conjTransposeGammaGai["Gai"] * (*GammaGai)["Gbj"];
    out.set<Tensor<complex> *>("PPHHCoulombIntegrals", Vabij);
  }

  if (Vaijb) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaijb->get_name() << std::endl;
    (*Vaijb)["aijb"] = conjTransposeGammaGai["Gaj"] * (*GammaGia)["Gib"];
    out.set<Tensor<complex> *>("PHHPCoulombIntegrals", Vaijb);
  }

  if (Vijab) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijab->get_name() << std::endl;
    (*Vijab)["ijab"] = conjTransposeGammaGia["Gia"] * (*GammaGia)["Gjb"];
    out.set<Tensor<complex> *>("HHPPCoulombIntegrals", Vijab);
  }

  if (Vaibj) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"] = conjTransposeGammaGab["Gab"] * (*GammaGij)["Gij"];
    out.set<Tensor<complex> *>("PHPHCoulombIntegrals", Vaibj);
  }

  if (Vijkl) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"] = conjTransposeGammaGij["Gik"] * (*GammaGij)["Gjl"];
    out.set<Tensor<complex> *>("HHHHCoulombIntegrals", Vijkl);
  }

  if (Vijka) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"] = conjTransposeGammaGij["Gik"] * (*GammaGia)["Gja"];
    out.set<Tensor<complex> *>("HHHPCoulombIntegrals", Vijka);
  }

  if (Vaijk) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaijk->get_name() << std::endl;
    (*Vaijk)["aijk"] = conjTransposeGammaGai["Gaj"] * (*GammaGij)["Gik"];
    out.set<Tensor<complex> *>("PHHHCoulombIntegrals", Vaijk);
  }

  if (Vabcd) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"] = conjTransposeGammaGab["Gac"] * (*GammaGab)["Gbd"];
    out.set<Tensor<complex> *>("PPPPCoulombIntegrals", Vabcd);
  }

  if (Vaibc) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibc->get_name() << std::endl;
    (*Vaibc)["aibc"] = conjTransposeGammaGab["Gab"] * (*GammaGia)["Gic"];
    out.set<Tensor<complex> *>("PHPPCoulombIntegrals", Vaibc);
  }

  // Force integrals to be real
  // --------------------------------------------------------
  if (in.get<bool>("forceReal")) {
    LOG(1, "CoulombIntegrals") << "Forcing integrals to be real" << std::endl;

    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vabij)["abij"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vijab)["ijab"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vaijb)["aijb"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vaibj)["aibj"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vijkl)["ijkl"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vijka)["ijka"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vaijk)["aijk"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vabcd)["abcd"]);
    CTF::Transform<complex>(std::function<void(complex &)>(
        [](complex &v) { v.imag(0.0); }))((*Vaibc)["aibc"]);
  }

  // Create the rest of integrals from the already given ones
  // --------------------------------------------------------

  if (Vabic) {
    // vvov = h%v * vovv
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vabic->get_name()
                               << " using " << Vaibc->get_name() << std::endl;

    Vabic->sum(1.0, *Vaibc, "ciba", 0.0, "abic", fConj);
    if (antisymmetrize) {
      // vvvo = h * vovv
      Vabic->sum(-1.0, *Vaibc, "ciab", 1.0, "abic", fConj);
    }
    out.set<Tensor<complex> *>("PPHPCoulombIntegrals", Vabic);
  }
  if (Vabci) {
    // vvvo = h * vovv
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vabci->get_name()
                               << " using " << Vaibc->get_name() << std::endl;

    Vabci->sum(1.0, *Vaibc, "ciab", 0.0, "abci", fConj);
    if (antisymmetrize) {
      // vvov = h%v * vovv
      Vabci->sum(-1.0, *Vaibc, "ciba", 1.0, "abci", fConj);
    }
    out.set<Tensor<complex> *>("PPPHCoulombIntegrals", Vabci);
  }
  if (Vijak) {
    // oovo = h * vooo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Vijak->get_name()
                               << " using " << Vaijk->get_name() << std::endl;

    Vijak->sum(1.0, *Vaijk, "akij", 0.0, "ijak", fConj);
    if (antisymmetrize) {
      // ooov = e * ooov
      Vijak->sum(-1.0, *Vijka, "ijka", 1.0, "ijak");
    }
    out.set<Tensor<complex> *>("HHPHCoulombIntegrals", Vijak);
  }
  if (Viajk) {
    // ovoo = h * ooov
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viajk->get_name()
                               << " using " << Vijka->get_name() << std::endl;

    Viajk->sum(1.0, *Vijka, "jkia", 0.0, "iajk", fConj);
    if (antisymmetrize) {
      // ovoo = h * ooov
      Viajk->sum(-1.0, *Vijka, "kjia", 1.0, "iajk", fConj);
    }
    out.set<Tensor<complex> *>("HPHHCoulombIntegrals", Viajk);
  }
  if (Viajb) {
    // ovov = v * vovo
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viajb->get_name()
                               << " using " << Vaibj->get_name() << std::endl;

    Viajb->sum(1.0, *Vaibj, "aibj", 0.0, "iajb");
    if (antisymmetrize) {
      // ovvo = h * voov
      Viajb->sum(-1.0, *Vaijb, "bjia", 1.0, "iajb", fConj);
    }
    out.set<Tensor<complex> *>("HPHPCoulombIntegrals", Viajb);
  }
  if (Viabj) {
    // ovvo = h * voov
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viabj->get_name()
                               << " using " << Vaijb->get_name() << std::endl;

    Viabj->sum(1.0, *Vaijb, "bjia", 0.0, "iabj", fConj);
    if (antisymmetrize) {
      // ovov = v * vovo
      Viabj->sum(-1.0, *Vaibj, "aibj", 1.0, "iabj");
    }
    out.set<Tensor<complex> *>("HPPHCoulombIntegrals", Viabj);
  }
  if (Viabc) {
    // ovvv = v * vovv
    LOG(1, "CoulombIntegrals") << "Evaluating " << Viabc->get_name()
                               << " using " << Vaibc->get_name() << std::endl;

    Viabc->sum(1.0, *Vaibc, "aicb", 0.0, "iabc");
    if (antisymmetrize) {
      // ovvv = v * vovv
      Viabc->sum(-1.0, *Vaibc, "aibc", 1.0, "iabc");
    }
    out.set<Tensor<complex> *>("HPPPCoulombIntegrals", Viabc);
  }

  if (antisymmetrize) {
    // Antisymmetrize integrals calculated directly from Gamma
    // IMPORTANT: This must be written after the creation of the integrals
    //            depending on them.
    if (Vijkl) (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
    if (Vijka) (*Vijka)["ijka"] -= (*Vijka)["jika"];
    if (Vabcd) (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
    if (Vabij) (*Vabij)["abij"] -= (*Vabij)["abji"];
    if (Vijab) (*Vijab)["ijab"] -= (*Vijab)["ijba"];
    if (Vaijk) (*Vaijk)["aijk"] -= (*Vaijk)["aikj"];
    if (Vaibc) (*Vaibc)["aibc"] -= (*Vaibc)["aicb"];

    // There is an inter-dependence of Vaijb and Vaibj for antisymmetrizing
    // so we define a temporary tensor, that is not antisymmetrized.
    Tensor<complex> originalVaijb(*Vaijb);

    if (Vaijb) (*Vaijb)["aijb"] -= (*Vaibj)["aibj"];
    if (Vaibj) (*Vaibj)["aibj"] -= originalVaijb["aijb"];
  }
}

void CoulombIntegralsFromVertex::dryCalculateComplexIntegrals() {
  DryTensor<complex> *Vaibj(
      out.present("PHPHCoulombIntegrals")
          ? new DryTensor<complex>(4, vovo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vabij(
      out.present("PPHHCoulombIntegrals")
          ? new DryTensor<complex>(4, vvoo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijab(
      out.present("HHPPCoulombIntegrals")
          ? new DryTensor<complex>(4, oovv.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vaijb(
      out.present("PHHPCoulombIntegrals")
          ? new DryTensor<complex>(4, voov.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijkl(
      out.present("HHHHCoulombIntegrals")
          ? new DryTensor<complex>(4, oooo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijka(
      out.present("HHHPCoulombIntegrals")
          ? new DryTensor<complex>(4, ooov.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vaijk(
      out.present("PHHHCoulombIntegrals")
          ? new DryTensor<complex>(4, vooo.data(), syms.data())
          : nullptr);

  if (Vaibj) { out.set<DryTensor<complex> *>("PHPHCoulombIntegrals", Vaibj); }
  if (Vabij) { out.set<DryTensor<complex> *>("PPHHCoulombIntegrals", Vabij); }
  if (Vijab) { out.set<DryTensor<complex> *>("HHPPCoulombIntegrals", Vijab); }
  if (Vaijb) { out.set<DryTensor<complex> *>("PHHPCoulombIntegrals", Vaijb); }
  if (Vijkl) { out.set<DryTensor<complex> *>("HHHHCoulombIntegrals", Vijkl); }
  if (Vijka) { out.set<DryTensor<complex> *>("HHHPCoulombIntegrals", Vijka); }
  if (Vaijk) { out.set<DryTensor<complex> *>("PHHHCoulombIntegrals", Vaijk); }
}
