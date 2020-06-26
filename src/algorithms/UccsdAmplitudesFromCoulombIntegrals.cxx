#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <unistd.h>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(UccsdAmplitudesFromCoulombIntegrals);

UccsdAmplitudesFromCoulombIntegrals::UccsdAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {}


UccsdAmplitudesFromCoulombIntegrals::~UccsdAmplitudesFromCoulombIntegrals() {
}

void UccsdAmplitudesFromCoulombIntegrals::run() {

  // Set this for cluster singles doubles algorithm
  setIntegerArgument("antisymmetrize", 1);
  setIntegerArgument("unrestricted", 1);

  checkArgumentsOrDie(
    { "intermediates" // stanton intermediates
    , "antisymmetrize", "unrestricted"
    // Solver
    , "mixer", "maxIterations", "amplitudesConvergence", "energyConvergence"
    , "mixingRatio", "maxResidua"
    // Fock matrix
    , "HPFockMatrix" , "PPFockMatrix" , "HHFockMatrix"
    , "HoleEigenEnergies" , "ParticleEigenEnergies"
    // integrals
    , "HHHHCoulombIntegrals" , "PPPPCoulombIntegrals" , "HHHPCoulombIntegrals"
    , "HHPPCoulombIntegrals" , "HPHHCoulombIntegrals" , "HPHPCoulombIntegrals"
    , "HPPPCoulombIntegrals" , "PPHHCoulombIntegrals" , "PPHPCoulombIntegrals"
    , "HPPHCoulombIntegrals" , "PHPPCoulombIntegrals" , "HHPHCoulombIntegrals"
    , "PPPHCoulombIntegrals" , "PHPHCoulombIntegrals" , "PHHPCoulombIntegrals"
    , "UccsdDoublesAmplitudes", "UccsdSinglesAmplitudes", "UccsdEnergy"
    }
  );

  usingIntermediates = getIntegerArgument("intermediates", 1) == 1;
  if (! usingIntermediates ) {
    LOG(0, getAbbreviation())
      << "Not using Stanton et al. intermediates, the code will be much slower."
      << std::endl;
  }
  LOG(0, getAbbreviation())
    << "stanton intermediates: " << usingIntermediates << std::endl;

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
  CTF::Tensor<F> *Fab, *Fij, *Fia;

  if (  isArgumentGiven("HPFockMatrix")
     && isArgumentGiven("HHFockMatrix")
     && isArgumentGiven("PPFockMatrix")
     )
  {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals "
                                << "since you provided FockMatrices\n";
    }
    Fia = getTensorArgument<F, CTF::Tensor<F> >("HPFockMatrix");
    Fab = getTensorArgument<F, CTF::Tensor<F> >("PPFockMatrix");
    Fij = getTensorArgument<F, CTF::Tensor<F> >("HHFockMatrix");
  } else {
    auto epsi = getTensorArgument<double>("HoleEigenEnergies")
       , epsa = getTensorArgument<double>("ParticleEigenEnergies")
       ;
    Fia = nullptr;
    const int Nv(epsa->lens[0])
            , No(epsi->lens[0])
            , vv[] = {Nv, Nv}
            , oo[] = {No, No}
            , syms[] = {NS, NS}
            ;
    // TODO: replace with unique_ptr and give away the hamiltonian one
    // such pointer
    Fab = new CTF::Tensor<F>(2, vv, syms, *Cc4s::world, "Fab");
    Fij = new CTF::Tensor<F>(2, oo, syms, *Cc4s::world, "Fij");
    CTF::Transform<double, F>([](double eps, F &f) { f = eps; })( (*epsi)["i"]
                                                                , (*Fij)["ii"]
                                                                );
    CTF::Transform<double, F>([](double eps, F &f) { f = eps; })( (*epsa)["a"]
                                                                , (*Fab)["aa"]
                                                                );
  }

  // Get couloumb integrals
  auto Vijkl(getTensorArgument<F>("HHHHCoulombIntegrals"))
     , Vabcd(getTensorArgument<F>("PPPPCoulombIntegrals"))
     , Vijka(getTensorArgument<F>("HHHPCoulombIntegrals"))
     , Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"))
     , Viajk(getTensorArgument<F>("HPHHCoulombIntegrals"))
     , Viajb(getTensorArgument<F>("HPHPCoulombIntegrals"))
     , Viabc(getTensorArgument<F>("HPPPCoulombIntegrals"))
     , Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"))
     , Vabic(getTensorArgument<F>("PPHPCoulombIntegrals"))
     , Viabj(getTensorArgument<F>("HPPHCoulombIntegrals"))
     , Vaibc(getTensorArgument<F>("PHPPCoulombIntegrals"))
     , Vijak(getTensorArgument<F>("HHPHCoulombIntegrals"))
     , Vabci(getTensorArgument<F>("PPPHCoulombIntegrals"))
     , Vaibj(getTensorArgument<F>("PHPHCoulombIntegrals"))
     , Vaijb(getTensorArgument<F>("PHHPCoulombIntegrals"))
     ;

  // Read the amplitudes Tai and Tabij
  auto Tai(amplitudes->get(0))
     , Tabij(amplitudes->get(1))
     ;

  auto residuum(NEW(FockVector<F>, *amplitudes));
  auto Rai(residuum->get(0))
     , Rabij(residuum->get(1))
     ;

  *residuum *= 0.0;
  Rai->set_name("Rai");
  Rabij->set_name("Rabij");

  if (iterationStep == 0){
    LOG(1, getAbbreviation())
      << "Set initial Rabij amplitudes to Vijab" << std::endl;
    (*Rabij)["abij"] = (*Vijab)["ijab"];
    return residuum;
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H
   // set fock matrix
   .setFij(Fij).setFab(Fab).setFia(Fia)
   // set coulomb integrals
   .setVabcd(Vabcd).setViajb(Viajb).setVijab(Vijab).setVijkl(Vijkl)
   .setVijka(Vijka).setViabc(Viabc).setViajk(Viajk).setVabic(Vabic)
   .setVaibc(Vaibc).setVaibj(Vaibj).setViabj(Viabj).setVijak(Vijak)
   .setVaijb(Vaijb).setVabci(Vabci).setVabij(Vabij)
   // set current t-amplitudes
   .setTai(Tai.get()).setTabij(Tabij.get())
   // set a general dressing, since we don't want any terms to get dropped
   .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL)
   // use stanton intermediates?
   .useStantonIntermediatesUCCSD(usingIntermediates)
   ;

  /* T1 equations:
   * =============
   *
   * The singles amplitude equations are simply taking the
   * Wai part of the \bar H and setting it to zero
   */
  auto Wai = H.getAI();
  (*Rai)["ai"]  = (*Wai)["ai"];
  //These are the residum equations, we have to substract them from Wai
  (*Rai)["bi"] += ( - 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Rai)["bi"] += ( + 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];

  /*
   * T2 equations:
   * =============
   *
   * The doubles amplitude equations are simply taking the
   * Wabij part of the \bar H and setting it to zero
   */
  auto Wabij = H.getABIJ();
  (*Rabij)["abij"]  = (*Wabij)["abij"];
  //These are the residum equations, substract them from Wabij
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];

  return residuum;

}
