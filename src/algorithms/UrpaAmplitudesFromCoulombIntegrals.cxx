#include <algorithms/UrpaAmplitudesFromCoulombIntegrals.hpp>
#include <equations/SimilarityTransformedHamiltonian.hpp>
#include <unistd.h>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>

using namespace sisi4s;

PTR(FockVector<sisi4s::complex>)
UrpaAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<complex>) &amplitudes) {
  if (iterationStep == 0) {
    LOG(1, getAbbreviation())
        << "WARNING: Using complex version of Urpa" << std::endl;
    LOG(1, getAbbreviation())
        << "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UrpaAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UrpaAmplitudesFromCoulombIntegrals::getResiduumTemplate(
    const int iterationStep,
    const PTR(const FockVector<F>) &amplitudes) {
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));

  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));

  // Get couloumb integrals
  auto Vijab(in.get<Tensor<F> *>("HHPPCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia;

  if (isArgumentGiven("HPFockMatrix") && isArgumentGiven("HHFockMatrix")
      && isArgumentGiven("PPFockMatrix")) {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    Fia = in.get<Tensor<F> *>("HPFockMatrix");
    Fab = in.get<Tensor<F> *>("PPFockMatrix");
    Fij = in.get<Tensor<F> *>("HHFockMatrix");
  } else {
    Fia = NULL;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");

  if (iterationStep == 0) {
    LOG(1, getAbbreviation())
        << "Set initial Rabij amplitudes to Vijab" << std::endl;
    (*Rabij)["abij"] = (*Vijab)["ijab"];
    return residuum;
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H.setFij(Fij)
      .setFab(Fab)
      .setFia(Fia)
      .setVijab(Vijab)
      .setTai(Tai.get())
      .setTabij(Tabij.get())
      .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL);

  // T1 equations:
  //
  // The singles amplitude equations are simply taking the
  // Wai part of the \bar H and setting it to zero
  //
  auto Wai = H.getAI_RPA();
  (*Rai)["ai"] += (*Wai)["ai"];
  // These are the residum equations, we have to substract them from Wai
  (*Rai)["bi"] += (-1.0) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Rai)["bi"] += (+1.0) * (*Fij)["ki"] * (*Tai)["bk"];

  // T2 equations:
  //
  // The doubles amplitude equations are simply taking the
  // Wabij part of the \bar H and setting it to zero
  //
  auto Wabij = H.getABIJ_RPA();
  (*Rabij)["abij"] += (*Wabij)["abij"];
  // These are the residum equations, substract them from Wabij
  (*Rabij)["cdij"] += (+1.0) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  (*Rabij)["cdij"] += (-1.0) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  (*Rabij)["cdij"] += (+1.0) * (*Fab)["de"] * (*Tabij)["ecij"];
  (*Rabij)["cdij"] += (-1.0) * (*Fab)["ce"] * (*Tabij)["edij"];

  return residuum;
}
