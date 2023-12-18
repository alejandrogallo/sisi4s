#include <algorithms/CcsdtEquationOfMotionDavidson.hpp>
#include <equations/SimilarityTransformedHamiltonian.hpp>
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

DEFSPEC(
    CcsdtEquationOfMotionDavidson,
    SPEC_IN(
        {"amplitudesConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},
        {"energyConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},
        {"preconditionerRandomSigma", SPEC_VALUE_DEF("TODO: DOC", double, 0.1)},
        {"complexVersion", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"eigenStates", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"intermediates", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"maxIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 32)},
        {"minIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"preconditionerRandom", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
        {"refreshOnMaxBasisSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
        {"oneBodyRdmRange", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
        {"refreshIterations", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
        {"DoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"ParticleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PHPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"SinglesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT());

IMPLEMENT_ALGORITHM(CcsdtEquationOfMotionDavidson) {

  if (in.get<int64_t>("complexVersion", 1) == 1) {
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
  bool preconditionerRandom(in.get<int64_t>("preconditionerRandom", 0) == 1);
  double preconditionerRandomSigma(
      in.get<double>("preconditionerRandomSigma", 0.1));
  bool refreshOnMaxBasisSize(in.get<int64_t>("refreshOnMaxBasisSize", 0) == 1);
  std::vector<int> oneBodyRdmIndices(
      RangeParser(in.get<std::string>("oneBodyRdmRange", "")).getRange());
  int eigenStates(in.get<int64_t>("eigenStates", 1));
  const double energyConvergence(in.get<double>("energyConvergence", 1e-6)),
      amplitudesConvergence(in.get<double>("amplitudesConvergence", 1e-6));
  bool intermediates(in.get<int64_t>("intermediates", 1));
  unsigned int maxIterations(in.get<int64_t>("maxIterations", 32));
  unsigned int minIterations(in.get<int64_t>("minIterations", 1));
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));
  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));
  std::vector<int> refreshIterations(
      RangeParser(in.get<std::string>("refreshIterations", "")).getRange());
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int maxBasisSize(
      in.get<int64_t>("maxBasisSize",
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
  Tensor<double> *pVijkl(in.get<Tensor<double> *>("HHHHCoulombIntegrals"));
  Tensor<F> cVijkl(pVijkl->order,
                   pVijkl->lens,
                   pVijkl->sym,
                   *Sisi4s::world,
                   pVijkl->get_name());
  Tensor<F> *Vijkl(&cVijkl);
  toComplexTensor(*pVijkl, *Vijkl);

  Tensor<double> *pVabcd(in.get<Tensor<double> *>("PPPPCoulombIntegrals"));
  Tensor<F> cVabcd(pVabcd->order,
                   pVabcd->lens,
                   pVabcd->sym,
                   *Sisi4s::world,
                   pVabcd->get_name());
  Tensor<F> *Vabcd(&cVabcd);
  toComplexTensor(*pVabcd, *Vabcd);

  Tensor<double> *pVijka(in.get<Tensor<double> *>("HHHPCoulombIntegrals"));
  Tensor<F> cVijka(pVijka->order,
                   pVijka->lens,
                   pVijka->sym,
                   *Sisi4s::world,
                   pVijka->get_name());
  Tensor<F> *Vijka(&cVijka);
  toComplexTensor(*pVijka, *Vijka);

  Tensor<double> *pVijab(in.get<Tensor<double> *>("HHPPCoulombIntegrals"));
  Tensor<F> cVijab(pVijab->order,
                   pVijab->lens,
                   pVijab->sym,
                   *Sisi4s::world,
                   pVijab->get_name());
  Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  Tensor<double> *pViajk(in.get<Tensor<double> *>("HPHHCoulombIntegrals"));
  Tensor<F> cViajk(pViajk->order,
                   pViajk->lens,
                   pViajk->sym,
                   *Sisi4s::world,
                   pViajk->get_name());
  Tensor<F> *Viajk(&cViajk);
  toComplexTensor(*pViajk, *Viajk);

  Tensor<double> *pViajb(in.get<Tensor<double> *>("HPHPCoulombIntegrals"));
  Tensor<F> cViajb(pViajb->order,
                   pViajb->lens,
                   pViajb->sym,
                   *Sisi4s::world,
                   pViajb->get_name());
  Tensor<F> *Viajb(&cViajb);
  toComplexTensor(*pViajb, *Viajb);

  Tensor<double> *pViabc(in.get<Tensor<double> *>("HPPPCoulombIntegrals"));
  Tensor<F> cViabc(pViabc->order,
                   pViabc->lens,
                   pViabc->sym,
                   *Sisi4s::world,
                   pViabc->get_name());
  Tensor<F> *Viabc(&cViabc);
  toComplexTensor(*pViabc, *Viabc);

  Tensor<double> *pVabic(in.get<Tensor<double> *>("PPHPCoulombIntegrals"));
  Tensor<F> cVabic(pVabic->order,
                   pVabic->lens,
                   pVabic->sym,
                   *Sisi4s::world,
                   pVabic->get_name());
  Tensor<F> *Vabic(&cVabic);
  toComplexTensor(*pVabic, *Vabic);

  Tensor<double> *pVabci(in.get<Tensor<double> *>("PPPHCoulombIntegrals"));
  Tensor<F> cVabci(pVabci->order,
                   pVabci->lens,
                   pVabci->sym,
                   *Sisi4s::world,
                   pVabci->get_name());
  Tensor<F> *Vabci(&cVabci);
  toComplexTensor(*pVabci, *Vabci);

  Tensor<double> *pVaibc(in.get<Tensor<double> *>("PHPPCoulombIntegrals"));
  Tensor<F> cVaibc(pVaibc->order,
                   pVaibc->lens,
                   pVaibc->sym,
                   *Sisi4s::world,
                   pVaibc->get_name());
  Tensor<F> *Vaibc(&cVaibc);
  toComplexTensor(*pVaibc, *Vaibc);

  Tensor<double> *pVaibj(in.get<Tensor<double> *>("PHPHCoulombIntegrals"));
  Tensor<F> cVaibj(pVaibj->order,
                   pVaibj->lens,
                   pVaibj->sym,
                   *Sisi4s::world,
                   pVaibj->get_name());
  Tensor<F> *Vaibj(&cVaibj);
  toComplexTensor(*pVaibj, *Vaibj);

  Tensor<double> *pViabj(in.get<Tensor<double> *>("HPPHCoulombIntegrals"));
  Tensor<F> cViabj(pViabj->order,
                   pViabj->lens,
                   pViabj->sym,
                   *Sisi4s::world,
                   pViabj->get_name());
  Tensor<F> *Viabj(&cViabj);
  toComplexTensor(*pViabj, *Viabj);

  Tensor<double> *pVijak(in.get<Tensor<double> *>("HHPHCoulombIntegrals"));
  Tensor<F> cVijak(pVijak->order,
                   pVijak->lens,
                   pVijak->sym,
                   *Sisi4s::world,
                   pVijak->get_name());
  Tensor<F> *Vijak(&cVijak);
  toComplexTensor(*pVijak, *Vijak);

  Tensor<double> *pVaijb(in.get<Tensor<double> *>("PHHPCoulombIntegrals"));
  Tensor<F> cVaijb(pVaijb->order,
                   pVaijb->lens,
                   pVaijb->sym,
                   *Sisi4s::world,
                   pVaijb->get_name());
  Tensor<F> *Vaijb(&cVaijb);
  toComplexTensor(*pVaijb, *Vaijb);

  // Tensor<double> *Vabij(
  // in.get<Tensor<double>*>("PPHHCoulombIntegrals"));

  // HF terms
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia(new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (isArgumentGiven("HPFockMatrix") && isArgumentGiven("HHFockMatrix")
      && isArgumentGiven("PPFockMatrix")) {
    LOG(0, "CcsdtEomDavid") << "Using non-canonical orbitals" << std::endl;

    Tensor<double> *realFia(in.get<Tensor<double> *>("HPFockMatrix"));
    Tensor<double> *realFab(in.get<Tensor<double> *>("PPFockMatrix"));
    Tensor<double> *realFij(in.get<Tensor<double> *>("HHFockMatrix"));
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
  toComplexTensor((*in.get<Tensor<double> *>("SinglesAmplitudes")), Tai);
  toComplexTensor((*in.get<Tensor<double> *>("DoublesAmplitudes")), Tabij);

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
