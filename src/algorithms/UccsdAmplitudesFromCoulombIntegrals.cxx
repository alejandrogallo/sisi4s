#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
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

// TODO remove F and FSPEC
#define F double
#define FSPEC double
DEFSPEC(UccsdAmplitudesFromCoulombIntegrals,
        SPEC_IN(UCCSD_SPEC_IN),
        SPEC_OUT(UCCSD_SPEC_OUT));
#undef F

IMPLEMENT_ALGORITHM(UccsdAmplitudesFromCoulombIntegrals) {

  // TODO: set this as an attribute in the ClusterSinglesDoublesAlgorithm
  // Set this for cluster singles doubles algorithm
  in.set<bool>("unrestricted", true);

  usingIntermediates = in.get<bool>("intermediates");
  onlyPpl = in.get<bool>("OnlyPPL");
  if (!usingIntermediates) {
    LOG(0, getAbbreviation()) << "Not using Stanton et al. intermediates, the "
                                 "code will be much slower."
                              << std::endl;
  }
  LOG(0, getAbbreviation())
      << "stanton intermediates: " << usingIntermediates << std::endl;

  ClusterSinglesDoublesAlgorithm::run();
}

PTR(FockVector<sisi4s::complex>)
UccsdAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<complex>) &amplitudes) {
  if (iterationStep == 0) {
    LOG(1, getAbbreviation())
        << "WARNING: Using complex version of Uccsd" << std::endl;
    LOG(1, getAbbreviation())
        << "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdAmplitudesFromCoulombIntegrals::getResiduumTemplate(
    const int iterationStep,
    const PTR(const FockVector<F>) &amplitudes) {
  const auto Vpqrs = in.get<CoulombIntegrals<F> *>("CoulombIntegrals");
  Tensor<F> *Fab, *Fij, *Fia;

  if (in.present("HPFockMatrix") && in.present("HHFockMatrix")
      && in.present("PPFockMatrix")) {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals "
                                << "since you provided FockMatrices\n";
    }
    Fia = in.get<Tensor<F> *>("HPFockMatrix");
    Fab = in.get<Tensor<F> *>("PPFockMatrix");
    Fij = in.get<Tensor<F> *>("HHFockMatrix");
  } else {
    auto epsi = in.get<Tensor<double> *>("HoleEigenEnergies"),
         epsa = in.get<Tensor<double> *>("ParticleEigenEnergies");
    Fia = nullptr;
    const int Nv(epsa->lens[0]), No(epsi->lens[0]),
        vv[] = {Nv, Nv}, oo[] = {No, No}, syms[] = {NS, NS};
    // TODO: replace with unique_ptr and give away the hamiltonian one
    // such pointer
    Fab = new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab");
    Fij = new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij");
    CTF::Transform<double, F>([](double eps, F &f) { f = eps; })((*epsi)["i"],
                                                                 (*Fij)["ii"]);
    CTF::Transform<double, F>([](double eps, F &f) { f = eps; })((*epsa)["a"],
                                                                 (*Fab)["aa"]);
  }

  const auto Vabij = Vpqrs->pphh();

  // Read the amplitudes Tai and Tabij
  auto Tai(amplitudes->get(0)), Tabij(amplitudes->get(1));

  auto residuum(NEW(FockVector<F>, *amplitudes));
  auto Rai(residuum->get(0)), Rabij(residuum->get(1));

  *residuum *= 0.0;
  Rai->set_name("Rai");
  Rabij->set_name("Rabij");

  if (iterationStep == 0) {
    if (onlyPpl) {
      const auto Vabcd = Vpqrs->pppp();
      LOG(1, "Performing only Ppl contraction") << std::endl;
      (*Rabij)["cdij"] += (+0.5) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
      (*Rabij)["cdij"] +=
          (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
      return residuum;
    } else {
      LOG(1, getAbbreviation())
          << "Set initial Rabij amplitudes to Vijab" << std::endl;
      (*Rabij)["abij"] = (*Vabij)["abij"];
      return residuum;
    }
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0], Vpqrs);

  H
      // set fock matrix
      .setFij(Fij)
      .setFab(Fab)
      .setFia(Fia)
      // set coulomb integrals
      // set current t-amplitudes
      .setTai(Tai.get())
      .setTabij(Tabij.get())
      // set a general dressing, since we don't want any terms to get dropped
      .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL)
      // use stanton intermediates?
      .withStantonIntermediatesUCCSD(usingIntermediates);

  // T1 equations:
  //
  // The singles amplitude equations are simply taking the
  // Wai part of the \bar H and setting it to zero
  //

  auto Wai = H.getAI();
  (*Rai)["ai"] = (*Wai)["ai"];
  // These are the residum equations, we have to substract them from Wai
  // NOTE: Notice that we're substracting only the diagonal of Fab and Fij
  (*Rai)["bi"] += (-1.0) * (*Fab)["bb"] * (*Tai)["bi"];
  (*Rai)["bi"] += (+1.0) * (*Fij)["ii"] * (*Tai)["bi"];

  /*
   * T2 equations:
   * =============
   *
   * The doubles amplitude equations are simply taking the
   * Wabij part of the \bar H and setting it to zero
   */
  auto Wabij = H.getABIJ();
  (*Rabij)["abij"] = (*Wabij)["abij"];
  // These are the residum equations, substract them from Wabij
  (*Rabij)["cdij"] += (+1.0) * (*Fij)["ii"] * (*Tabij)["cdij"];
  (*Rabij)["cdij"] += (-1.0) * (*Fij)["jj"] * (*Tabij)["cdji"];
  (*Rabij)["cdij"] += (+1.0) * (*Fab)["dd"] * (*Tabij)["dcij"];
  (*Rabij)["cdij"] += (-1.0) * (*Fab)["cc"] * (*Tabij)["cdij"];

  return residuum;
}
