#include <algorithms/UnrestrictedEquationOfMotionSinglesFromRpa.hpp>
#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <algorithms/CcsdPreconditioner.hpp>
#include <algorithms/OneBodyReducedDensityMatrix.hpp>

#include <math/EigenSystemDavidson.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/TensorIo.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(UnrestrictedEquationOfMotionSinglesFromRpa);

UnrestrictedEquationOfMotionSinglesFromRpa::UnrestrictedEquationOfMotionSinglesFromRpa(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
UnrestrictedEquationOfMotionSinglesFromRpa::~UnrestrictedEquationOfMotionSinglesFromRpa() {}


void UnrestrictedEquationOfMotionSinglesFromRpa::run() {

  if (getIntegerArgument("complexVersion", 1) == 1) {
    LOG(0, "RPAEomDavid") << "Using complex code" << std::endl;
    UnrestrictedEquationOfMotionSinglesFromRpa::run<complex>();
  } else {
    LOG(0, "RPAEomDavid") << "Using real code" << std::endl;
    UnrestrictedEquationOfMotionSinglesFromRpa::run<double>();
  }

}

template <typename F>
void UnrestrictedEquationOfMotionSinglesFromRpa::run() {

  // Arguments
  bool preconditionerRandom(
    getIntegerArgument("preconditionerRandom", 0) == 1
  );
  double preconditionerRandomSigma(getRealArgument(
    "preconditionerRandomSigma", 0.1
  ));
  bool refreshOnMaxBasisSize(
    getIntegerArgument("refreshOnMaxBasisSize", 0) == 1
  );
  std::vector<int> oneBodyRdmIndices(
    RangeParser(getTextArgument("oneBodyRdmRange", "")).getRange()
  );
  int eigenStates(getIntegerArgument("eigenstates", 1));
  double ediff(getRealArgument("ediff", 1e-4));
  bool intermediates(getIntegerArgument("intermediates", 1));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));
  std::vector<int> eigenvectorsIndices(
    RangeParser(getTextArgument("printEigenvectorsRange", "")).getRange()
  );
  CTF::Tensor<double> *epsi(
    getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies")
  );
  CTF::Tensor<double> *epsa(
    getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies")
  );
  std::vector<int> refreshIterations(
    RangeParser(getTextArgument("refreshIterations", "")).getRange()
  );
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int  maxBasisSize(getIntegerArgument(
    "maxBasisSize", No*Nv + (No*(No - 1)/2 ) * (Nv * (Nv - 1)/2)
  ));

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int vo[] = {Nv,No};
  int oo[] = {No, No};
  int vvoo[] = {Nv,Nv,No,No};

  // Logging arguments
  LOG(0, "RPAEomDavid") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "RPAEomDavid") << "ediff " << ediff << std::endl;
  LOG(0, "RPAEomDavid") << eigenStates << " eigen states" << std::endl;
  LOG(0, "RPAEomDavid") << "No: " << No << std::endl;
  LOG(0, "RPAEomDavid") << "Nv: " << Nv << std::endl;
  LOG(0, "RPAEomDavid") << "maxBasisSize: " << maxBasisSize << std::endl;


  // Get copy of couloumb integrals

  CTF::Tensor<double> *pVijab(
    getTensorArgument<double, CTF::Tensor<double> >("HHPPCoulombIntegrals")
  );
  CTF::Tensor<F> cVijab(
    pVijab->order, pVijab->lens, pVijab->sym, *Cc4s::world,
    pVijab->get_name()
  );
  CTF::Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  CTF::Tensor<double> *pViajb(
    getTensorArgument<double, CTF::Tensor<double> >("HPHPCoulombIntegrals")
  );
  CTF::Tensor<F> cViajb(
    pViajb->order, pViajb->lens, pViajb->sym, *Cc4s::world,
    pViajb->get_name()
  );
  CTF::Tensor<F> *Viajb(&cViajb);
  toComplexTensor(*pViajb, *Viajb);

  CTF::Tensor<double> *pVaibj(
    getTensorArgument<double, CTF::Tensor<double> >("PHPHCoulombIntegrals")
  );
  CTF::Tensor<F> cVaibj(
    pVaibj->order, pVaibj->lens, pVaibj->sym, *Cc4s::world,
    pVaibj->get_name()
  );
  CTF::Tensor<F> *Vaibj(&cVaibj);
  toComplexTensor(*pVaibj, *Vaibj);

  CTF::Tensor<double> *pViabj(
    getTensorArgument<double, CTF::Tensor<double> >("HPPHCoulombIntegrals")
  );
  CTF::Tensor<F> cViabj(
    pViabj->order, pViabj->lens, pViabj->sym, *Cc4s::world,
    pViabj->get_name()
  );
  CTF::Tensor<F> *Viabj(&cViabj);
  toComplexTensor(*pViabj, *Viabj);

  CTF::Tensor<double> *pVaijb(
    getTensorArgument<double, CTF::Tensor<double> >("PHHPCoulombIntegrals")
  );
  CTF::Tensor<F> cVaijb(
    pVaijb->order, pVaijb->lens, pVaijb->sym, *Cc4s::world,
    pVaijb->get_name()
  );
  CTF::Tensor<F> *Vaijb(&cVaijb);
  toComplexTensor(*pVaijb, *Vaijb);

  // HF terms
  CTF::Tensor<F> *Fab(new CTF::Tensor<F>(2, vv, syms2, *Cc4s::world, "Fab"));
  CTF::Tensor<F> *Fij(new CTF::Tensor<F>(2, oo, syms2, *Cc4s::world, "Fij"));
  CTF::Tensor<F> *Fia(new CTF::Tensor<F>(2, ov, syms2, *Cc4s::world, "Fia"));

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    LOG(0, "RPAEomDavid") << "Using non-canonical orbitals" << std::endl;

    CTF::Tensor<double> *realFia(
      getTensorArgument<double, CTF::Tensor<double> >("HPFockMatrix")
    );
    CTF::Tensor<double> *realFab(
      getTensorArgument<double, CTF::Tensor<double> >("PPFockMatrix")
    );
    CTF::Tensor<double> *realFij(
      getTensorArgument<double, CTF::Tensor<double> >("HHFockMatrix")
    );
    toComplexTensor(*realFij, *Fij);
    toComplexTensor(*realFab, *Fab);
    toComplexTensor(*realFia, *Fia);
  } else {
    LOG(0, "RPAEomDavid") << "Using canonical orbitals" << std::endl;
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

  CTF::Tensor<F> Tai(2, vo, syms2, *Cc4s::world, "Tai");
  CTF::Tensor<F> Tabij(4, vvoo, syms4, *Cc4s::world, "Tabij");
  toComplexTensor(
    (*getTensorArgument<double, CTF::Tensor<double> >("SinglesAmplitudes")),
    Tai
  );
  toComplexTensor(
    (*getTensorArgument<double, CTF::Tensor<double> >("DoublesAmplitudes")),
    Tabij
  );

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H.setFij(Fij).setFab(Fab).setFia(Fia)
    .setViajb(Viajb)
    .setVijab(Vijab)
    .setVaibj(Vaibj)
    .setViabj(Viabj)
    .setVaijb(Vaijb)
    .setTai(&Tai).setTabij(&Tabij)
    .setRightApplyIntermediates(intermediates)
    .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::RPA);

  struct RpaH {
  public:
    SimilarityTransformedHamiltonian<F>* h;
    SFockVector<F> rightApply(SFockVector<F> &V) {
      return h->rightApplyHirata_RPA(V);
    }
  } rpaH;
  rpaH.h = &H;

  // INITIALIZE SIMILARITY PRECONDITIONER
  CcsdPreconditioner<F> P;
  P
    .setTai(&Tai).setTabij(&Tabij)
    .setFij(Fij).setFab(Fab)
    // Set coulomb integrals
    .setViajb(Viajb).setVijab(Vijab)
    // Set random information
    .setRandom(preconditionerRandom).setRandomSigma(preconditionerRandomSigma)
  ;

  EigenSystemDavidsonMono<
    RpaH,
    CcsdPreconditioner<F>,
    SFockVector<F>
  > eigenSystem(
    &rpaH,
    eigenStates,
    &P,
    ediff,
    maxBasisSize,
    maxIterations,
    minIterations
  );
  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "RPAEomDavid") <<
      "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();

  if (eigenvectorsIndices.size() > 0) {

    for (auto &index: eigenvectorsIndices) {
      LOG(1, "RPAEomDavid") << "Writing out eigenvector " << index << std::endl;
      auto eigenState(eigenSystem.getRightEigenVectors()[index-1]);
      TensorIo::writeText<F>(
        "Rai-" + std::to_string(index) + ".tensor",
        *eigenState.get(0),
        "ij", "", " "
      );
    }
  }

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  NEW_FILE("EomRPAEnergies.dat") << "";
  for (auto &ev: eigenValues) {
    eigenCounter++;
    LOG(0, "RPAEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("EomRPAEnergies.dat") << eigenCounter <<
      " " << ev.real() << " " << ev.imag() << std::endl;
  }

}
