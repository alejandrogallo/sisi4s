#include <algorithms/CcsdtEquationOfMotionDavidson.hpp>
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
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/SharedPointer.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdtEquationOfMotionDavidson);

CcsdtEquationOfMotionDavidson::CcsdtEquationOfMotionDavidson(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}
CcsdtEquationOfMotionDavidson::~CcsdtEquationOfMotionDavidson() {}

void CcsdtEquationOfMotionDavidson::run() {

  if (getIntegerArgument("complexVersion", 1) == 1) {
    LOG(0, "CcsdtEomDavid") << "Using complex code" << std::endl;
    CcsdtEquationOfMotionDavidson::run<complex>();
  } else {
    LOG(0, "CcsdtEomDavid") << "Using real code" << std::endl;
    CcsdtEquationOfMotionDavidson::run<double>();
  }
}

template <typename F>
void CcsdtEquationOfMotionDavidson::run() {

  // Arguments
  bool preconditionerRandom(getIntegerArgument("preconditionerRandom", 0) == 1);
  double preconditionerRandomSigma(
      getRealArgument("preconditionerRandomSigma", 0.1));
  bool refreshOnMaxBasisSize(getIntegerArgument("refreshOnMaxBasisSize", 0)
                             == 1);
  std::vector<int> oneBodyRdmIndices(
      RangeParser(getTextArgument("oneBodyRdmRange", "")).getRange());
  int eigenStates(getIntegerArgument("eigenStates", 1));
  const double energyConvergence(getRealArgument("energyConvergence", 1e-6)),
      amplitudesConvergence(getRealArgument("amplitudesConvergence", 1e-6));
  bool intermediates(getIntegerArgument("intermediates", 1));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));
  Tensor<double> *epsi(
      getTensorArgument<double, Tensor<double>>("HoleEigenEnergies"));
  Tensor<double> *epsa(
      getTensorArgument<double, Tensor<double>>("ParticleEigenEnergies"));
  std::vector<int> refreshIterations(
      RangeParser(getTextArgument("refreshIterations", "")).getRange());
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int maxBasisSize(
      getIntegerArgument("maxBasisSize",
                         No * Nv + (No * (No - 1) / 2) * (Nv * (Nv - 1) / 2)));

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int vo[] = {Nv, No};
  int oo[] = {No, No};
  int vvoo[] = {Nv, Nv, No, No};

  // Logging arguments
  LOG(0, "CcsdtEomDavid") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "CcsdtEomDavid") << "energyConvergence " << energyConvergence
                          << std::endl;
  LOG(0, "CcsdtEomDavid") << eigenStates << " eigen states" << std::endl;
  LOG(0, "CcsdtEomDavid") << "No: " << No << std::endl;
  LOG(0, "CcsdtEomDavid") << "Nv: " << Nv << std::endl;
  LOG(0, "CcsdtEomDavid") << "maxBasisSize: " << maxBasisSize << std::endl;

  // Get copy of couloumb integrals
  Tensor<double> *pVijkl(
      getTensorArgument<double, Tensor<double>>("HHHHCoulombIntegrals"));
  Tensor<F> cVijkl(pVijkl->order,
                   pVijkl->lens,
                   pVijkl->sym,
                   *Sisi4s::world,
                   pVijkl->get_name());
  Tensor<F> *Vijkl(&cVijkl);
  toComplexTensor(*pVijkl, *Vijkl);

  Tensor<double> *pVabcd(
      getTensorArgument<double, Tensor<double>>("PPPPCoulombIntegrals"));
  Tensor<F> cVabcd(pVabcd->order,
                   pVabcd->lens,
                   pVabcd->sym,
                   *Sisi4s::world,
                   pVabcd->get_name());
  Tensor<F> *Vabcd(&cVabcd);
  toComplexTensor(*pVabcd, *Vabcd);

  Tensor<double> *pVijka(
      getTensorArgument<double, Tensor<double>>("HHHPCoulombIntegrals"));
  Tensor<F> cVijka(pVijka->order,
                   pVijka->lens,
                   pVijka->sym,
                   *Sisi4s::world,
                   pVijka->get_name());
  Tensor<F> *Vijka(&cVijka);
  toComplexTensor(*pVijka, *Vijka);

  Tensor<double> *pVijab(
      getTensorArgument<double, Tensor<double>>("HHPPCoulombIntegrals"));
  Tensor<F> cVijab(pVijab->order,
                   pVijab->lens,
                   pVijab->sym,
                   *Sisi4s::world,
                   pVijab->get_name());
  Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  Tensor<double> *pViajk(
      getTensorArgument<double, Tensor<double>>("HPHHCoulombIntegrals"));
  Tensor<F> cViajk(pViajk->order,
                   pViajk->lens,
                   pViajk->sym,
                   *Sisi4s::world,
                   pViajk->get_name());
  Tensor<F> *Viajk(&cViajk);
  toComplexTensor(*pViajk, *Viajk);

  Tensor<double> *pViajb(
      getTensorArgument<double, Tensor<double>>("HPHPCoulombIntegrals"));
  Tensor<F> cViajb(pViajb->order,
                   pViajb->lens,
                   pViajb->sym,
                   *Sisi4s::world,
                   pViajb->get_name());
  Tensor<F> *Viajb(&cViajb);
  toComplexTensor(*pViajb, *Viajb);

  Tensor<double> *pViabc(
      getTensorArgument<double, Tensor<double>>("HPPPCoulombIntegrals"));
  Tensor<F> cViabc(pViabc->order,
                   pViabc->lens,
                   pViabc->sym,
                   *Sisi4s::world,
                   pViabc->get_name());
  Tensor<F> *Viabc(&cViabc);
  toComplexTensor(*pViabc, *Viabc);

  Tensor<double> *pVabic(
      getTensorArgument<double, Tensor<double>>("PPHPCoulombIntegrals"));
  Tensor<F> cVabic(pVabic->order,
                   pVabic->lens,
                   pVabic->sym,
                   *Sisi4s::world,
                   pVabic->get_name());
  Tensor<F> *Vabic(&cVabic);
  toComplexTensor(*pVabic, *Vabic);

  Tensor<double> *pVabci(
      getTensorArgument<double, Tensor<double>>("PPPHCoulombIntegrals"));
  Tensor<F> cVabci(pVabci->order,
                   pVabci->lens,
                   pVabci->sym,
                   *Sisi4s::world,
                   pVabci->get_name());
  Tensor<F> *Vabci(&cVabci);
  toComplexTensor(*pVabci, *Vabci);

  Tensor<double> *pVaibc(
      getTensorArgument<double, Tensor<double>>("PHPPCoulombIntegrals"));
  Tensor<F> cVaibc(pVaibc->order,
                   pVaibc->lens,
                   pVaibc->sym,
                   *Sisi4s::world,
                   pVaibc->get_name());
  Tensor<F> *Vaibc(&cVaibc);
  toComplexTensor(*pVaibc, *Vaibc);

  Tensor<double> *pVaibj(
      getTensorArgument<double, Tensor<double>>("PHPHCoulombIntegrals"));
  Tensor<F> cVaibj(pVaibj->order,
                   pVaibj->lens,
                   pVaibj->sym,
                   *Sisi4s::world,
                   pVaibj->get_name());
  Tensor<F> *Vaibj(&cVaibj);
  toComplexTensor(*pVaibj, *Vaibj);

  Tensor<double> *pViabj(
      getTensorArgument<double, Tensor<double>>("HPPHCoulombIntegrals"));
  Tensor<F> cViabj(pViabj->order,
                   pViabj->lens,
                   pViabj->sym,
                   *Sisi4s::world,
                   pViabj->get_name());
  Tensor<F> *Viabj(&cViabj);
  toComplexTensor(*pViabj, *Viabj);

  Tensor<double> *pVijak(
      getTensorArgument<double, Tensor<double>>("HHPHCoulombIntegrals"));
  Tensor<F> cVijak(pVijak->order,
                   pVijak->lens,
                   pVijak->sym,
                   *Sisi4s::world,
                   pVijak->get_name());
  Tensor<F> *Vijak(&cVijak);
  toComplexTensor(*pVijak, *Vijak);

  Tensor<double> *pVaijb(
      getTensorArgument<double, Tensor<double>>("PHHPCoulombIntegrals"));
  Tensor<F> cVaijb(pVaijb->order,
                   pVaijb->lens,
                   pVaijb->sym,
                   *Sisi4s::world,
                   pVaijb->get_name());
  Tensor<F> *Vaijb(&cVaijb);
  toComplexTensor(*pVaijb, *Vaijb);

  // Tensor<double> *Vabij(
  // getTensorArgument<double, Tensor<double>>("PPHHCoulombIntegrals"));

  // HF terms
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia(new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (isArgumentGiven("HPFockMatrix") && isArgumentGiven("HHFockMatrix")
      && isArgumentGiven("PPFockMatrix")) {
    LOG(0, "CcsdtEomDavid") << "Using non-canonical orbitals" << std::endl;

    Tensor<double> *realFia(
        getTensorArgument<double, Tensor<double>>("HPFockMatrix"));
    Tensor<double> *realFab(
        getTensorArgument<double, Tensor<double>>("PPFockMatrix"));
    Tensor<double> *realFij(
        getTensorArgument<double, Tensor<double>>("HHFockMatrix"));
    toComplexTensor(*realFij, *Fij);
    toComplexTensor(*realFab, *Fab);
    toComplexTensor(*realFia, *Fia);
  } else {
    LOG(0, "CcsdtEomDavid") << "Using canonical orbitals" << std::endl;
    Fia = NULL;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  Tensor<F> Tai(2, vo, syms2, *Sisi4s::world, "Tai");
  Tensor<F> Tabij(4, vvoo, syms4, *Sisi4s::world, "Tabij");
  toComplexTensor(
      (*getTensorArgument<double, Tensor<double>>("SinglesAmplitudes")),
      Tai);
  toComplexTensor(
      (*getTensorArgument<double, Tensor<double>>("DoublesAmplitudes")),
      Tabij);

  enum SimilarityTransformedHamiltonian<F>::Dressing dressing;
  if (!isArgumentGiven("TriplesAmplitudes")) {
    dressing = SimilarityTransformedHamiltonian<F>::Dressing::CCSD;
  } else {
    dressing = SimilarityTransformedHamiltonian<F>::Dressing::CCSDT;
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H.setFij(Fij)
      .setFab(Fab)
      .setFia(Fia)
      .setVabcd(Vabcd)
      .setViajb(Viajb)
      .setVijab(Vijab)
      .setVijkl(Vijkl)
      .setVijka(Vijka)
      .setViabc(Viabc)
      .setViajk(Viajk)
      .setVabic(Vabic)
      .setVaibc(Vaibc)
      .setVaibj(Vaibj)
      .setViabj(Viabj)
      .setVijak(Vijak)
      .setVaijb(Vaijb)
      .setVabci(Vabci)
      .setTai(&Tai)
      .setTabij(&Tabij)
      .with_right_apply_intermediates(intermediates)
      .setDressing(dressing);

  CcsdPreconditioner<F>
      P(Tai, Tabij, *Fij, *Fab, *Vabcd, *Viajb, *Vijab, *Vijkl);
  P.preconditionerRandom = preconditionerRandom;
  P.preconditionerRandomSigma = preconditionerRandomSigma;

  EigenSystemDavidsonMono<SimilarityTransformedHamiltonian<F>,
                          CcsdPreconditioner<F>,
                          SDTFockVector<F>>
      eigenSystem(&H,
                  eigenStates,
                  &P,
                  amplitudesConvergence,
                  energyConvergence,
                  maxBasisSize,
                  maxIterations,
                  minIterations);
  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "CcsdtEomDavid")
        << "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  NEW_FILE("EomCcsdtEnergies.dat") << "";
  for (auto &ev : eigenValues) {
    eigenCounter++;
    LOG(0, "CcsdtEomDavid")
        << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("EomCcsdtEnergies.dat")
        << eigenCounter << " " << ev.real() << " " << ev.imag() << std::endl;
  }
}
