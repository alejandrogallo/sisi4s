#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <unistd.h>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(UccsdAmplitudesFromCoulombIntegrals);

UccsdAmplitudesFromCoulombIntegrals::UccsdAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {}


UccsdAmplitudesFromCoulombIntegrals::~UccsdAmplitudesFromCoulombIntegrals() {
}

void UccsdAmplitudesFromCoulombIntegrals::run() {
  usingIntermediates = (bool) getIntegerArgument("intermediates", 1);
  if (! usingIntermediates ) {
    LOG(0, getAbbreviation()) <<
      "Not using intermediates, the code will be much slower."
    << std::endl;
  }
  ClusterSinglesDoublesAlgorithm::run();
}

PTR(FockVector<cc4s::complex>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<complex>) &amplitudes
) {
  if (iterationStep == 0){
    LOG(1, getAbbreviation()) <<
       "WARNING: Using complex version of Uccsd" << std::endl;
    LOG(1, getAbbreviation()) <<
       "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<double>) &amplitudes
) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdAmplitudesFromCoulombIntegrals::getResiduumTemplate(
  const int iterationStep, const PTR(const FockVector<F>) &amplitudes
) {
  CTF::Tensor<double> *epsi(
    getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies")
  );

  CTF::Tensor<double> *epsa(
    getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies")
  );

  // Get couloumb integrals
  auto Vijkl(getTensorArgument<F, CTF::Tensor<F> >("HHHHCoulombIntegrals"));
  auto Vabcd(getTensorArgument<F, CTF::Tensor<F> >("PPPPCoulombIntegrals"));
  auto Vijka(getTensorArgument<F, CTF::Tensor<F> >("HHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<F, CTF::Tensor<F> >("HHPPCoulombIntegrals"));
  auto Viajk(getTensorArgument<F, CTF::Tensor<F> >("HPHHCoulombIntegrals"));
  auto Viajb(getTensorArgument<F, CTF::Tensor<F> >("HPHPCoulombIntegrals"));
  auto Viabc(getTensorArgument<F, CTF::Tensor<F> >("HPPPCoulombIntegrals"));
  auto Vabij(getTensorArgument<F, CTF::Tensor<F> >("PPHHCoulombIntegrals"));
  auto Vabic(getTensorArgument<F, CTF::Tensor<F> >("PPHPCoulombIntegrals"));
  auto Viabj(getTensorArgument<F, CTF::Tensor<F> >("HPPHCoulombIntegrals"));
  auto Vaibc(getTensorArgument<F, CTF::Tensor<F> >("PHPPCoulombIntegrals"));
  auto Vijak(getTensorArgument<F, CTF::Tensor<F> >("HHPHCoulombIntegrals"));
  auto Vabci(getTensorArgument<F, CTF::Tensor<F> >("PPPHCoulombIntegrals"));
  auto Vaibj(getTensorArgument<F, CTF::Tensor<F> >("PHPHCoulombIntegrals"));
  auto Vaijb(getTensorArgument<F, CTF::Tensor<F> >("PHHPCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  CTF::Tensor<F> *Fab(
    new CTF::Tensor<F>(2, vv, syms, *Cc4s::world, "Fab")
  );
  CTF::Tensor<F> *Fij(
    new CTF::Tensor<F>(2, oo, syms, *Cc4s::world, "Fij")
  );
  CTF::Tensor<F> *Fia;

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    if (iterationStep == 0){
    LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    Fia = getTensorArgument<F, CTF::Tensor<F> >("HPFockMatrix");
    Fab = getTensorArgument<F, CTF::Tensor<F> >("PPFockMatrix");
    Fij = getTensorArgument<F, CTF::Tensor<F> >("HHFockMatrix");
  } else {
    Fia = NULL;
    CTF::Transform<double, F>(
      std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }
      )
    ) (
      (*epsi)["i"], (*Fij)["ii"]
    );
    CTF::Transform<double, F>(
      std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }
      )
    ) (
      (*epsa)["a"], (*Fab)["aa"]
    );
  }

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  //amplitudes->get(0)
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

  if (iterationStep == 0){
    LOG(1, getAbbreviation())
      << "Set initial Rabij amplitudes to Vijab"
    << std::endl;
    (*Rabij)["abij"] = (*Vijab)["ijab"];
    return residuum;
  }

  SimilarityTransformedHamiltonian<F> H(
    Fij, Fab, Fia,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci, Vabij,
    usingIntermediates,
    SimilarityTransformedHamiltonian<F>::Dressing::GENERAL
  );
  H.setTai(Tai.get());
  H.setTabij(Tabij.get());

  // T1 equations:
  //
  // The singles amplitude equations are simply taking the
  // Wai part of the \bar H and setting it to zero
  //
  auto Wai = H.getAI();
  (*Rai)["ai"] += (*Wai)["ai"];
  //These are the residum equations, we have to substract them from Wai
  (*Rai)["bi"] += ( - 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Rai)["bi"] += ( + 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];

  // T2 equations:
  //
  // The doubles amplitude equations are simply taking the
  // Wabij part of the \bar H and setting it to zero
  //
  auto Wabij = H.getABIJ();
  (*Rabij)["abij"] += (*Wabij)["abij"];
  //These are the residum equations, substract them from Wabij
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];

  return residuum;

}
