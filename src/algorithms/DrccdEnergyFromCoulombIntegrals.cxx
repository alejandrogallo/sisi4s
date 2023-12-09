#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

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
  auto Vabij(in.get<Tensor<F> *>("PPHHCoulombIntegrals"));
  auto Vaijb(in.get<Tensor<F> *>("PHHPCoulombIntegrals"));
  auto Vijab(in.get<Tensor<F> *>("HHPPCoulombIntegrals"));

  // Check for spin polarization
  double spins(in.get<int64_t>("unrestricted", 0) ? 1.0 : 2.0);
  // get amplitude parts
  auto Tabij(amplitudes->get(1));

  // construct residuum
  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= F(0);
  auto Rabij(residuum->get(1));

  int linearized(in.get<int64_t>("linearized", 0));
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
      if (in.get<int64_t>("adjacentPairsExchange", 0)) {
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
