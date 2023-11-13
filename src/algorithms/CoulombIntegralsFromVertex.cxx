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

ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromVertex);

CoulombIntegralsFromVertex::CoulombIntegralsFromVertex(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

CoulombIntegralsFromVertex::~CoulombIntegralsFromVertex() {}

void CoulombIntegralsFromVertex::run() {
  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<double> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<double> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

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
  if (isArgumentGiven("PHPHCoulombIntegrals")) {
    list.push_back("PHPHCoulombIntegrals");
  }
  if (isArgumentGiven("PPHHCoulombIntegrals")) {
    list.push_back("PPHHCoulombIntegrals");
  }
  if (isArgumentGiven("HHHHCoulombIntegrals")) {
    list.push_back("HHHHCoulombIntegrals");
  }
  if (isArgumentGiven("HHHPCoulombIntegrals")) {
    list.push_back("HHHPCoulombIntegrals");
  }
  if (isArgumentGiven("PPPPCoulombIntegrals")) {
    list.push_back("PPPPCoulombIntegrals");
  }
  if (isArgumentGiven("PPPHCoulombIntegrals")) {
    list.push_back("PPPHCoulombIntegrals");
  }
  if (isArgumentGiven("PHHHCoulombIntegrals")) {
    list.push_back("PHHHCoulombIntegrals");
  }
  if (isArgumentGiven("HHPPCoulombIntegrals")) {
    list.push_back("HHPPCoulombIntegrals");
  }
  if (isArgumentGiven("PHHPCoulombIntegrals")) {
    list.push_back("PHHPCoulombIntegrals");
  }
  if (isArgumentGiven("HPHHCoulombIntegrals")) {
    list.push_back("HPHHCoulombIntegrals");
  }
  if (isArgumentGiven("HPHPCoulombIntegrals")) {
    list.push_back("HPHPCoulombIntegrals");
  }
  if (isArgumentGiven("HPPPCoulombIntegrals")) {
    list.push_back("HPPPCoulombIntegrals");
  }
  if (isArgumentGiven("PPHPCoulombIntegrals")) {
    list.push_back("PPHPCoulombIntegrals");
  }
  if (isArgumentGiven("HPPHCoulombIntegrals")) {
    list.push_back("HPPHCoulombIntegrals");
  }
  if (isArgumentGiven("HHPHCoulombIntegrals")) {
    list.push_back("HHPHCoulombIntegrals");
  }
  if (isArgumentGiven("PHPPCoulombIntegrals")) {
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

  bool realIntegrals = !getIntegerArgument("complex", 0);
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
  DryTensor<complex> *GammaGqr(
      getTensorArgument<complex, DryTensor<complex>>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
      getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(
      getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));

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

  bool realIntegrals = !getIntegerArgument("complex", 0);
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

  const bool calculate_vpqrs = isArgumentGiven("PQRSCoulombIntegrals");

  int antisymmetrize(getIntegerArgument("antisymmetrize", 0));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals")
        << "Calculating antisymmetrized integrals" << std::endl;
  }
  Tensor<double> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vovo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaibj")
                            : nullptr);
  Tensor<double> *Vabij(isArgumentGiven("PPHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvoo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabij")
                            : nullptr);
  Tensor<double> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oooo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijkl")
                            : nullptr);
  Tensor<double> *Vijka(isArgumentGiven("HHHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ooov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijka")
                            : nullptr);
  Tensor<double> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvvv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabcd")
                            : nullptr);
  Tensor<double> *Vabci(isArgumentGiven("PPPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvvo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabci")
                            : nullptr);

  // Initialization of tensors created from already existing ones
  Tensor<double> *Vaijk(isArgumentGiven("PHHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vooo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaijk")
                            : nullptr);
  Tensor<double> *Vijab(isArgumentGiven("HHPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oovv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijab")
                            : nullptr);
  Tensor<double> *Vaijb(isArgumentGiven("PHHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 voov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vaijb")
                            : nullptr);
  Tensor<double> *Viajk(isArgumentGiven("HPHHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovoo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viajk")
                            : nullptr);
  Tensor<double> *Viajb(isArgumentGiven("HPHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viajb")
                            : nullptr);
  Tensor<double> *Viabc(isArgumentGiven("HPPPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovvv.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viabc")
                            : nullptr);
  Tensor<double> *Vabic(isArgumentGiven("PPHPCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 vvov.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vabic")
                            : nullptr);
  Tensor<double> *Viabj(isArgumentGiven("HPPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 ovvo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Viabj")
                            : nullptr);
  Tensor<double> *Vijak(isArgumentGiven("HHPHCoulombIntegrals")
                            ? new Tensor<double>(4,
                                                 oovo.data(),
                                                 syms.data(),
                                                 *Sisi4s::world,
                                                 "Vijak")
                            : nullptr);
  Tensor<double> *Vaibc(isArgumentGiven("PHPPCoulombIntegrals")
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
    Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
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
    allocatedTensorArgument("PQRSCoulombIntegrals", t);
  }

  // Compute the integrals Vabij Vaibj Vaijb Vijkl Vabcd
  if (Vaibj) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"] = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
    allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabij->get_name() << std::endl;
    (*Vabij)["abij"] = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijkl) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"] = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
    allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"] = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
    allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  }
  if (Vabcd) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"] = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
    allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vabci) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabci->get_name() << std::endl;
    (*Vabci)["abci"] = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
    allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);
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
    allocatedTensorArgument("HPPHCoulombIntegrals", Viabj);
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
    allocatedTensorArgument("HPHPCoulombIntegrals", Viajb);
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
    allocatedTensorArgument("HPPPCoulombIntegrals", Viabc);
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
    allocatedTensorArgument("HHPHCoulombIntegrals", Vijak);
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
    allocatedTensorArgument("HHPPCoulombIntegrals", Vijab);
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
    allocatedTensorArgument("PPHPCoulombIntegrals", Vabic);
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
    allocatedTensorArgument("PHHPCoulombIntegrals", Vaijb);
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
    allocatedTensorArgument("PHPPCoulombIntegrals", Vaibc);
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
    allocatedTensorArgument("PHHHCoulombIntegrals", Vaijk);
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
    allocatedTensorArgument("HPHHCoulombIntegrals", Viajk);
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
  DryTensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals")
                         ? new DryTensor<>(4, vovo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals")
                         ? new DryTensor<>(4, vvoo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijab(isArgumentGiven("HHPPCoulombIntegrals")
                         ? new DryTensor<>(4, oovv.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vaijb(isArgumentGiven("PHHPCoulombIntegrals")
                         ? new DryTensor<>(4, voov.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals")
                         ? new DryTensor<>(4, oooo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals")
                         ? new DryTensor<>(4, ooov.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vaijk(isArgumentGiven("PHHHCoulombIntegrals")
                         ? new DryTensor<>(4, vooo.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals")
                         ? new DryTensor<>(4, vvvv.data(), syms.data())
                         : nullptr);
  DryTensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals")
                         ? new DryTensor<>(4, vvvo.data(), syms.data())
                         : nullptr);

  if (Vaibj) { allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj); }
  if (Vabij) { allocatedTensorArgument("PPHHCoulombIntegrals", Vabij); }
  if (Vijab) { allocatedTensorArgument("HHPPCoulombIntegrals", Vijab); }
  if (Vaijb) { allocatedTensorArgument("PHHPCoulombIntegrals", Vaijb); }
  if (Vijkl) { allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl); }
  if (Vijka) { allocatedTensorArgument("HHHPCoulombIntegrals", Vijka); }
  if (Vaijk) { allocatedTensorArgument("PHHHCoulombIntegrals", Vaijk); }
  if (Vabcd) { allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd); }
  if (Vabci) { allocatedTensorArgument("PPPHCoulombIntegrals", Vabci); }
}

void CoulombIntegralsFromVertex::calculateComplexIntegrals() {
  int antisymmetrize(getIntegerArgument("antisymmetrize", 0));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals")
        << "Calculating antisymmetrized integrals" << std::endl;
  }

  Tensor<complex> *Vabij(isArgumentGiven("PPHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvoo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabij")
                             : nullptr);

  Tensor<complex> *Vijab(
      // TODO: HHPP is always conj(Permute(PPHH))
      isArgumentGiven("HHPPCoulombIntegrals")
          ? new Tensor<complex>(4,
                                oovv.data(),
                                syms.data(),
                                *Sisi4s::world,
                                "Vijab")
          : nullptr);

  Tensor<complex> *Vaijb(isArgumentGiven("PHHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   voov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaijb")
                             : nullptr);

  Tensor<complex> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vovo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaibj")
                             : nullptr);

  Tensor<complex> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   oooo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijkl")
                             : nullptr);

  Tensor<complex> *Vijka(isArgumentGiven("HHHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ooov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijka")
                             : nullptr);

  Tensor<complex> *Vaijk(
      // TODO: PHHH is always conj(Permute(HHHP))
      isArgumentGiven("PHHHCoulombIntegrals")
          ? new Tensor<complex>(4,
                                vooo.data(),
                                syms.data(),
                                *Sisi4s::world,
                                "Vaijk")
          : nullptr);

  Tensor<complex> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvvv.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabcd")
                             : nullptr);

  Tensor<complex> *Vaibc(isArgumentGiven("PHPPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vovv.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vaibc")
                             : nullptr);

  // Initialization of tensors created from already existing ones
  // In principle the integrals above do not constitute the minimal set of
  // integrals from which one can write the rest (which would be 7)

  Tensor<complex> *Vabic(isArgumentGiven("PPHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabic")
                             : nullptr);
  Tensor<complex> *Vabci(isArgumentGiven("PPPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   vvvo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vabci")
                             : nullptr);
  Tensor<complex> *Vijak(isArgumentGiven("HHPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   oovo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Vijak")
                             : nullptr);
  Tensor<complex> *Viajk(isArgumentGiven("HPHHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovoo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viajk")
                             : nullptr);
  Tensor<complex> *Viajb(isArgumentGiven("HPHPCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovov.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viajb")
                             : nullptr);
  Tensor<complex> *Viabj(isArgumentGiven("HPPHCoulombIntegrals")
                             ? new Tensor<complex>(4,
                                                   ovvo.data(),
                                                   syms.data(),
                                                   *Sisi4s::world,
                                                   "Viabj")
                             : nullptr);
  Tensor<complex> *Viabc(isArgumentGiven("HPPPCoulombIntegrals")
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
    allocatedTensorArgument<complex>("PPHHCoulombIntegrals", Vabij);
  }

  if (Vaijb) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaijb->get_name() << std::endl;
    (*Vaijb)["aijb"] = conjTransposeGammaGai["Gaj"] * (*GammaGia)["Gib"];
    allocatedTensorArgument<complex>("PHHPCoulombIntegrals", Vaijb);
  }

  if (Vijab) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijab->get_name() << std::endl;
    (*Vijab)["ijab"] = conjTransposeGammaGia["Gia"] * (*GammaGia)["Gjb"];
    allocatedTensorArgument<complex>("HHPPCoulombIntegrals", Vijab);
  }

  if (Vaibj) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"] = conjTransposeGammaGab["Gab"] * (*GammaGij)["Gij"];
    allocatedTensorArgument<complex>("PHPHCoulombIntegrals", Vaibj);
  }

  if (Vijkl) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"] = conjTransposeGammaGij["Gik"] * (*GammaGij)["Gjl"];
    allocatedTensorArgument<complex>("HHHHCoulombIntegrals", Vijkl);
  }

  if (Vijka) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"] = conjTransposeGammaGij["Gik"] * (*GammaGia)["Gja"];
    allocatedTensorArgument<complex>("HHHPCoulombIntegrals", Vijka);
  }

  if (Vaijk) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaijk->get_name() << std::endl;
    (*Vaijk)["aijk"] = conjTransposeGammaGai["Gaj"] * (*GammaGij)["Gik"];
    allocatedTensorArgument<complex>("PHHHCoulombIntegrals", Vaijk);
  }

  if (Vabcd) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"] = conjTransposeGammaGab["Gac"] * (*GammaGab)["Gbd"];
    allocatedTensorArgument<complex>("PPPPCoulombIntegrals", Vabcd);
  }

  if (Vaibc) {
    LOG(1, "CoulombIntegrals")
        << "Evaluating " << Vaibc->get_name() << std::endl;
    (*Vaibc)["aibc"] = conjTransposeGammaGab["Gab"] * (*GammaGia)["Gic"];
    allocatedTensorArgument<complex>("PHPPCoulombIntegrals", Vaibc);
  }

  // Force integrals to be real
  // --------------------------------------------------------
  if (getIntegerArgument("forceReal", 0) == 1) {
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
    allocatedTensorArgument<complex>("PPHPCoulombIntegrals", Vabic);
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
    allocatedTensorArgument<complex>("PPPHCoulombIntegrals", Vabci);
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
    allocatedTensorArgument<complex>("HHPHCoulombIntegrals", Vijak);
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
    allocatedTensorArgument<complex>("HPHHCoulombIntegrals", Viajk);
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
    allocatedTensorArgument<complex>("HPHPCoulombIntegrals", Viajb);
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
    allocatedTensorArgument<complex>("HPPHCoulombIntegrals", Viabj);
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
    allocatedTensorArgument<complex>("HPPPCoulombIntegrals", Viabc);
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
      isArgumentGiven("PHPHCoulombIntegrals")
          ? new DryTensor<complex>(4, vovo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vabij(
      isArgumentGiven("PPHHCoulombIntegrals")
          ? new DryTensor<complex>(4, vvoo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijab(
      isArgumentGiven("HHPPCoulombIntegrals")
          ? new DryTensor<complex>(4, oovv.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vaijb(
      isArgumentGiven("PHHPCoulombIntegrals")
          ? new DryTensor<complex>(4, voov.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijkl(
      isArgumentGiven("HHHHCoulombIntegrals")
          ? new DryTensor<complex>(4, oooo.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vijka(
      isArgumentGiven("HHHPCoulombIntegrals")
          ? new DryTensor<complex>(4, ooov.data(), syms.data())
          : nullptr);
  DryTensor<complex> *Vaijk(
      isArgumentGiven("PHHHCoulombIntegrals")
          ? new DryTensor<complex>(4, vooo.data(), syms.data())
          : nullptr);

  if (Vaibj) {
    allocatedTensorArgument<complex, DryTensor<complex>>("PHPHCoulombIntegrals",
                                                         Vaibj);
  }
  if (Vabij) {
    allocatedTensorArgument<complex, DryTensor<complex>>("PPHHCoulombIntegrals",
                                                         Vabij);
  }
  if (Vijab) {
    allocatedTensorArgument<complex, DryTensor<complex>>("HHPPCoulombIntegrals",
                                                         Vijab);
  }
  if (Vaijb) {
    allocatedTensorArgument<complex, DryTensor<complex>>("PHHPCoulombIntegrals",
                                                         Vaijb);
  }
  if (Vijkl) {
    allocatedTensorArgument<complex, DryTensor<complex>>("HHHHCoulombIntegrals",
                                                         Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument<complex, DryTensor<complex>>("HHHPCoulombIntegrals",
                                                         Vijka);
  }
  if (Vaijk) {
    allocatedTensorArgument<complex, DryTensor<complex>>("PHHHCoulombIntegrals",
                                                         Vaijk);
  }
}
