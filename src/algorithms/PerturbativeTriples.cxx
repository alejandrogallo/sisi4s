#include <algorithms/PerturbativeTriples.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriples);

PerturbativeTriples::PerturbativeTriples(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

PerturbativeTriples::~PerturbativeTriples() {}

void PerturbativeTriples::runInMemory() {
  Tensor<double> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<double> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<double> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<double> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<double> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<double> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<double> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  int vvvooo[] = {Nv, Nv, Nv, No, No, No};
  int syms[] = {NS, NS, NS, NS, NS, NS};
  Tensor<double> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");

  {
    Tensor<double> Zabcijk(6, vvvooo, syms, *Vabij->wrld, "Zabcijk");
    Zabcijk["abcijk"] = (*Vabci)["bcdk"] * (*Tabij)["adij"];
    Zabcijk["abcijk"] -= (*Vijka)["jklc"] * (*Tabij)["abil"];

    Tabcijk["abcijk"] = (*epsi)["i"];
    Tabcijk["abcijk"] += (*epsi)["j"];
    Tabcijk["abcijk"] += (*epsi)["k"];
    Tabcijk["abcijk"] -= (*epsa)["a"];
    Tabcijk["abcijk"] -= (*epsa)["b"];
    Tabcijk["abcijk"] -= (*epsa)["c"];
    Bivar_Function<> fDivide(&divide<double>);
    Tabcijk.contract(1.0,
                     Zabcijk,
                     "abcijk",
                     Tabcijk,
                     "abcijk",
                     0.0,
                     "abcijk",
                     fDivide);

    Zabcijk["abcijk"] = Tabcijk["abcijk"];
    Tabcijk["abcijk"] += Zabcijk["bacjik"];
    Tabcijk["abcijk"] += Zabcijk["acbikj"];
    Tabcijk["abcijk"] += Zabcijk["cbakji"];
    Tabcijk["abcijk"] += Zabcijk["cabkij"];
    Tabcijk["abcijk"] += Zabcijk["bcajki"];
  }

  Tensor<double> Zai(false, *Tai);
  Zai.set_name("Zai");
  Tensor<double> Zabij(false, *Tabij);
  Zabij.set_name("Zabij");

  Zai["ai"] = (+2.0) * Tabcijk["acdikl"] * (*Vabij)["cdkl"];
  Zai["ai"] += (-1.0) * Tabcijk["acdikl"] * (*Vabij)["dckl"];
  Zai["ai"] += (-2.0) * Tabcijk["acdlki"] * (*Vabij)["cdkl"];
  Zai["ai"] += (+1.0) * Tabcijk["acdlki"] * (*Vabij)["dckl"];

  Zabij["abij"] = (+2.0) * Tabcijk["acdijk"] * (*Vabci)["cdbk"];
  Zabij["abij"] += (-1.0) * Tabcijk["acdijk"] * (*Vabci)["dcbk"];
  Zabij["abij"] += (-1.0) * Tabcijk["acdkji"] * (*Vabci)["cdbk"];

  Zabij["abij"] += (-2.0) * Tabcijk["abcikl"] * (*Vijka)["kljc"];
  Zabij["abij"] += (+1.0) * Tabcijk["abcikl"] * (*Vijka)["lkjc"];
  Zabij["abij"] += (+1.0) * Tabcijk["abclki"] * (*Vijka)["kljc"];

  Scalar<> energy(*Sisi4s::world);
  double e, triplese;
  double ccsde(getRealArgument("CcsdEnergy"));

  energy[""] = (+2.0) * Zai["ai"] * (*Tai)["ai"];
  energy[""] += (+4.0) * Zabij["abij"] * (*Tabij)["abij"];
  energy[""] += (-2.0) * Zabij["abij"] * (*Tabij)["abji"];
  triplese = energy.get_val();
  e = triplese + ccsde;

  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << ccsde << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << triplese << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::runPiecuch() {
  Tensor<double> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<double> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<double> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<double> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<double> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<double> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<double> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  int vvvooo[] = {Nv, Nv, Nv, No, No, No};
  int syms[] = {NS, NS, NS, NS, NS, NS};
  Tensor<double> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");
  Tensor<double> Xabcijk(6, vvvooo, syms, *Vabij->wrld, "Xabcijk");

  Tabcijk["abcijk"] = (*Vabci)["bcek"] * (*Tabij)["aeij"];
  Tabcijk["abcijk"] -= (*Vijka)["jkmc"] * (*Tabij)["abim"];

  Xabcijk["abcijk"] = Tabcijk["abcijk"];
  Tabcijk["abcijk"] += Xabcijk["bacjik"];
  Tabcijk["abcijk"] += Xabcijk["acbikj"];
  Tabcijk["abcijk"] += Xabcijk["cbakji"];
  Tabcijk["abcijk"] += Xabcijk["cabkij"];
  Tabcijk["abcijk"] += Xabcijk["bcajki"];

  {
    Tensor<double> Zabcijk(6, vvvooo, syms, *Vabij->wrld, "Zabcijk");
    Zabcijk["abcijk"] = (1.0) * (*Tai)["ai"] * (*Vabij)["bcjk"];
    Zabcijk["abcijk"] += (1.0) * (*Tai)["bj"] * (*Vabij)["acik"];
    Zabcijk["abcijk"] += (1.0) * (*Tai)["ck"] * (*Vabij)["abij"];

    Xabcijk["abcijk"] = (4.0 / 3.0) * Zabcijk["abcijk"];
    Xabcijk["abcijk"] += (-2.0) * Zabcijk["acbijk"];
    Xabcijk["abcijk"] += (2.0 / 3.0) * Zabcijk["bcaijk"];

    Xabcijk["abcijk"] += (4.0 / 3.0) * Tabcijk["abcijk"];
    Xabcijk["abcijk"] += (-2.0) * Tabcijk["acbijk"];
    Xabcijk["abcijk"] += (2.0 / 3.0) * Tabcijk["bcaijk"];

    Zabcijk["abcijk"] = (*epsi)["i"];
    Zabcijk["abcijk"] += (*epsi)["j"];
    Zabcijk["abcijk"] += (*epsi)["k"];
    Zabcijk["abcijk"] -= (*epsa)["a"];
    Zabcijk["abcijk"] -= (*epsa)["b"];
    Zabcijk["abcijk"] -= (*epsa)["c"];
    Bivar_Function<> fDivide(&divide<double>);
    Tabcijk.contract(1.0,
                     Tabcijk,
                     "abcijk",
                     Zabcijk,
                     "abcijk",
                     0.0,
                     "abcijk",
                     fDivide);
  }

  Scalar<> energy(*Sisi4s::world);
  double e, triplese;
  double ccsde(getRealArgument("CcsdEnergy"));

  energy[""] = Xabcijk["abcijk"] * Tabcijk["abcijk"];
  triplese = energy.get_val();
  e = triplese + ccsde;

  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << ccsde << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << triplese << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::run() {
  Tensor<double> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<double> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<double> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<double> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<double> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<double> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<double> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  int vvvooo[] = {Nv, Nv, Nv, No, No, No};
  int syms[] = {NS, NS, NS, NS, NS, NS};
  Tensor<double> SVabcijk(6, vvvooo, syms, *Vabij->wrld, "SVabcijk");
  SVabcijk["abcijk"] = 0.5 * (*Tai)["ai"] * (*Vabij)["bcjk"];

  Tensor<double> DVabcijk(6, vvvooo, syms, *Vabij->wrld, "DVabcijk");
  DVabcijk["abcijk"] = (*Vabci)["bcdk"] * (*Tabij)["adij"];
  DVabcijk["abcijk"] -= (*Vijka)["jklc"] * (*Tabij)["abil"];

  Tensor<double> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");
  Tabcijk["abcijk"] = (+8.0) * DVabcijk["abcijk"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["acbijk"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["bacijk"];
  Tabcijk["abcijk"] += (+2.0) * DVabcijk["bcaijk"];
  Tabcijk["abcijk"] += (+2.0) * DVabcijk["cabijk"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["cbaijk"];

  Tabcijk["abcijk"] += (+8.0) * SVabcijk["abcijk"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["acbijk"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["bacijk"];
  Tabcijk["abcijk"] += (+2.0) * SVabcijk["bcaijk"];
  Tabcijk["abcijk"] += (+2.0) * SVabcijk["cabijk"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["cbaijk"];

  SVabcijk["abcijk"] = (*epsi)["i"];
  SVabcijk["abcijk"] += (*epsi)["j"];
  SVabcijk["abcijk"] += (*epsi)["k"];
  SVabcijk["abcijk"] -= (*epsa)["a"];
  SVabcijk["abcijk"] -= (*epsa)["b"];
  SVabcijk["abcijk"] -= (*epsa)["c"];
  Bivar_Function<> fDivide(&divide<double>);
  Tabcijk.contract(1.0,
                   Tabcijk,
                   "abcijk",
                   SVabcijk,
                   "abcijk",
                   0.0,
                   "abcijk",
                   fDivide);

  Scalar<> energy(*Sisi4s::world);
  energy[""] = DVabcijk["abcijk"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["bacjik"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["acbikj"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["cbakji"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["cabkij"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["bcajki"] * Tabcijk["abcijk"];

  double eTriples(energy.get_val());
  double eCcsd(getRealArgument("CcsdEnergy"));
  double e(eCcsd + eTriples);
  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << eTriples << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PPPHCoulombIntegrals");

  DryTensor<> *Tai(
      getTensorArgument<double, DryTensor<double>>("CcsdSinglesAmplitudes"));
  DryTensor<> *Tabij(
      getTensorArgument<double, DryTensor<double>>("CcsdDoublesAmplitudes"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
      getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(
      getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int vvvooo[] = {Nv, Nv, Nv, No, No, No};
  int syms[] = {NS, NS, NS, NS, NS, NS};
  DryTensor<> Tabcijk(6, vvvooo, syms, SOURCE_LOCATION);

  { DryTensor<> Zabcijk(6, vvvooo, syms, SOURCE_LOCATION); }

  DryTensor<> Zai(*Tai, SOURCE_LOCATION);
  DryTensor<> Zabij(*Tabij, SOURCE_LOCATION);

  DryScalar<> energy();
}

void PerturbativeTriples::dryRunPiecuch() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PPPHCoulombIntegrals");

  getTensorArgument<double, DryTensor<double>>("CcsdSinglesAmplitudes");
  getTensorArgument<double, DryTensor<double>>("CcsdDoublesAmplitudes");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
      getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(
      getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int vvvooo[] = {Nv, Nv, Nv, No, No, No};
  int syms[] = {NS, NS, NS, NS, NS, NS};
  DryTensor<> Tabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  DryTensor<> Yabcijk(6, vvvooo, syms, SOURCE_LOCATION);

  {
    DryTensor<> Zabcijk(6, vvvooo, syms, SOURCE_LOCATION);
    DryTensor<> Xabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  }

  DryScalar<> energy();
}
