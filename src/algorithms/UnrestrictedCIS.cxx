#include <algorithms/CcsdPreconditioner.hpp>
#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <algorithms/UnrestrictedCIS.hpp>

#include <Sisi4s.hpp>
#include <math/ComplexTensor.hpp>
#include <math/EigenSystemDavidson.hpp>
#include <math/FockVector.hpp>
#include <math/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <util/RangeParser.hpp>
#include <util/RangeParser.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>
#include <util/TensorIo.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace sisi4s;

enum Mode { CISD, CIS };

template <Mode m, typename F>
static void runGeneral(Algorithm *alg) {

  // Arguments
  bool preconditionerRandom(alg->getIntegerArgument("preconditionerRandom", 0)
                            == 1);
  double preconditionerRandomSigma(
      alg->getRealArgument("preconditionerRandomSigma", 0.1));
  bool refreshOnMaxBasisSize(alg->getIntegerArgument("refreshOnMaxBasisSize", 0)
                             == 1);
  std::vector<int> oneBodyRdmIndices(
      RangeParser(alg->getTextArgument("oneBodyRdmRange", "")).getRange());
  int eigenStates(alg->getIntegerArgument("eigenStates", 1));
  const double energyConvergence(
      alg->getRealArgument("energyConvergence", 1e-6)),
      amplitudesConvergence(
          alg->getRealArgument("amplitudesConvergence", 1e-6));
  bool intermediates(alg->getIntegerArgument("intermediates", 1));
  unsigned int maxIterations(alg->getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(alg->getIntegerArgument("minIterations", 1));
  std::vector<int> eigenvectorsIndices(
      RangeParser(alg->getTextArgument("printEigenvectorsRange", ""))
          .getRange());
  Tensor<double> *epsi(
      alg->getTensorArgument<double, Tensor<double>>("HoleEigenEnergies"));
  Tensor<double> *epsa(
      alg->getTensorArgument<double, Tensor<double>>("ParticleEigenEnergies"));
  std::vector<int> refreshIterations(
      RangeParser(alg->getTextArgument("refreshIterations", "")).getRange());
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int maxBasisSize(alg->getIntegerArgument(
      "maxBasisSize",
      No * Nv + (No * (No - 1) / 2) * (Nv * (Nv - 1) / 2)));

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int vo[] = {Nv, No};
  int oo[] = {No, No};
  int vvoo[] = {Nv, Nv, No, No};

  // Logging arguments
  LOG(0, "uCIS") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "uCIS") << "energyConvergence " << energyConvergence << std::endl;
  LOG(0, "uCIS") << eigenStates << " eigen states" << std::endl;
  LOG(0, "uCIS") << "No: " << No << std::endl;
  LOG(0, "uCIS") << "Nv: " << Nv << std::endl;
  LOG(0, "uCIS") << "maxBasisSize: " << maxBasisSize << std::endl;

  // Get copy of couloumb integrals

  Tensor<double> *pVijab(
      alg->getTensorArgument<double, Tensor<double>>("HHPPCoulombIntegrals"));
  Tensor<F> cVijab(pVijab->order,
                   pVijab->lens,
                   pVijab->sym,
                   *Sisi4s::world,
                   pVijab->get_name());
  Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  Tensor<double> *pViajb(
      alg->getTensorArgument<double, Tensor<double>>("HPHPCoulombIntegrals"));
  Tensor<F> cViajb(pViajb->order,
                   pViajb->lens,
                   pViajb->sym,
                   *Sisi4s::world,
                   pViajb->get_name());
  Tensor<F> *Viajb(&cViajb);
  toComplexTensor(*pViajb, *Viajb);

  // HF terms
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia(new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (alg->isArgumentGiven("HPFockMatrix")
      && alg->isArgumentGiven("HHFockMatrix")
      && alg->isArgumentGiven("PPFockMatrix")) {
    LOG(0, "uCIS") << "Using non-canonical orbitals" << std::endl;

    Tensor<double> *realFia(
        alg->getTensorArgument<double, Tensor<double>>("HPFockMatrix"));
    Tensor<double> *realFab(
        alg->getTensorArgument<double, Tensor<double>>("PPFockMatrix"));
    Tensor<double> *realFij(
        alg->getTensorArgument<double, Tensor<double>>("HHFockMatrix"));
    toComplexTensor(*realFij, *Fij);
    toComplexTensor(*realFab, *Fab);
    toComplexTensor(*realFia, *Fia);
  } else {
    LOG(0, "uCIS") << "Using canonical orbitals" << std::endl;
    Fia = NULL;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  if (m == CIS) {
    H.setFij(Fij)
        .setFab(Fab)
        .setFia(Fia)
        .setVijab(Vijab)
        .setViajb(Viajb)
        .withCIS(true)
        .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL);
  } else if (m == CISD) {
    throw "Not implemented CISD";
    H.setFij(Fij)
        .setFab(Fab)
        .setFia(Fia)
        .setVijab(Vijab)
        .setViajb(Viajb)
        .withCISD(true)
        .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL);
  }

  struct RpaH {
  public:
    SimilarityTransformedHamiltonian<F> *h;
    SFockVector<F> right_apply(SFockVector<F> &V) {
      return h->right_apply_CISD(V);
    }
  } rpaH;
  rpaH.h = &H;

  // INITIALIZE SIMILARITY PRECONDITIONER
  CISPreconditioner<F> P;
  P.setFij(Fij)
      .setFab(Fab)
      // Set coulomb integrals
      .setVijab(Vijab)
      // Set random information
      .setRandom(preconditionerRandom)
      .setRandomSigma(preconditionerRandomSigma);

  EigenSystemDavidsonMono<RpaH, CISPreconditioner<F>, SFockVector<F>>
      eigenSystem(&rpaH,
                  eigenStates,
                  &P,
                  amplitudesConvergence,
                  energyConvergence,
                  maxBasisSize,
                  maxIterations,
                  minIterations);
  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "uCIS") << "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();

  if (eigenvectorsIndices.size() > 0) {

    for (auto &index : eigenvectorsIndices) {
      LOG(1, "uCIS") << "Writing out eigenvector " << index << std::endl;
      auto eigenState(eigenSystem.getRightEigenVectors()[index - 1]);
      TensorIo::writeText<F>("Rai-" + std::to_string(index) + ".tensor",
                             *eigenState.get(0),
                             "ij",
                             "",
                             " ");
    }
  }

  std::vector<sisi4s::complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  NEW_FILE("EomRPAEnergies.dat") << "";
  for (auto &ev : eigenValues) {
    eigenCounter++;
    LOG(0, "uCIS") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("CISEnergies.dat")
        << eigenCounter << " " << ev.real() << " " << ev.imag() << std::endl;
  }
}

IMPLEMENT_ALGORITHM(UnrestrictedCIS) {
  if (getIntegerArgument("complexVersion", 1) == 1) {
    LOG(0, "uCIS") << "Using complex code" << std::endl;
    runGeneral<CIS, complex>(this);
  } else {
    LOG(0, "uCIS") << "Using real code" << std::endl;
    runGeneral<CIS, double>(this);
  }
}
