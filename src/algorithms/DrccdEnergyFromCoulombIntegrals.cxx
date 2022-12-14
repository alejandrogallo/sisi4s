#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombIntegrals);

DrccdEnergyFromCoulombIntegrals::DrccdEnergyFromCoulombIntegrals(
    std::vector<Argument> const &argumentList)
    : ClusterSinglesDoublesAlgorithm(argumentList) {}

DrccdEnergyFromCoulombIntegrals::~DrccdEnergyFromCoulombIntegrals() {}

PTR(FockVector<double>) DrccdEnergyFromCoulombIntegrals::getResiduum(
    const int iteration,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduum<double>(iteration, amplitudes);
}

PTR(FockVector<sisi4s::complex>) DrccdEnergyFromCoulombIntegrals::getResiduum(
    const int iteration,
    const PTR(const FockVector<sisi4s::complex>) &amplitudes) {
  return getResiduum<sisi4s::complex>(iteration, amplitudes);
}

template <typename F>
PTR(FockVector<F>) DrccdEnergyFromCoulombIntegrals::getResiduum(
    const int iteration,
    const PTR(const FockVector<F>) &amplitudes) {
  // read all required integrals
  auto Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
  auto Vaijb(getTensorArgument<F>("PHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));

  // Check for spin polarization
  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  // get amplitude parts
  auto Tabij(amplitudes->get(1));

  // construct residuum
  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= F(0);
  auto Rabij(residuum->get(1));

  int linearized(getIntegerArgument("linearized", 0));
  if (linearized) {
    LOG(1, getCapitalizedAbbreviation())
        << "Solving linearized T2 Amplitude Equations" << std::endl;
  } else {
    LOG(1, getCapitalizedAbbreviation())
        << "Solving T2 Amplitude Equations" << std::endl;
  }

  if (iteration > 0 || isArgumentGiven("startingDoublesAmplitudes")) {
    // for the remaining iterations compute the drCCD residuum
    (*Rabij)["abij"] += (*Vabij)["abij"];
    (*Rabij)["abij"] += spins * (*Vaijb)["akic"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += spins * (*Vaijb)["bkjc"] * (*Tabij)["acik"];
    if (!linearized) {
      Tensor<F> Wijab(false, *Vijab);
      Wijab["ijab"] = spins * (*Vijab)["ijab"];
      if (getIntegerArgument("adjacentPairsExchange", 0)) {
        Wijab["ijab"] -= (*Vijab)["jiab"];
      }
      // Construct intermediates
      Tensor<F> Calid(false, *Vaijb);
      Calid["alid"] = spins * Wijab["klcd"] * (*Tabij)["acik"];
      (*Rabij)["abij"] += Calid["alid"] * (*Tabij)["dblj"];
    }
  } else {
    // no amplitudes given: start with MP2 amplitudes
    (*Rabij)["abij"] += (*Vabij)["abij"];
  }

  return residuum;
}
