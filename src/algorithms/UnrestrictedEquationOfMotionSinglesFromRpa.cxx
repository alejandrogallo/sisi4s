#include <algorithms/UnrestrictedEquationOfMotionSinglesFromRpa.hpp>
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

IMPLEMENT_EMPTY_DRYRUN(UnrestrictedEquationOfMotionSinglesFromRpa) {}

DEFSPEC(
    UnrestrictedEquationOfMotionSinglesFromRpa,
    SPEC_IN(
        {"amplitudesConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},
        {"energyConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},
        {"preconditionerRandomSigma", SPEC_VALUE_DEF("TODO: DOC", double, 0.1)},
        {"complexVersion", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"eigenstates", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"intermediates", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"maxIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 32)},
        {"minIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"preconditionerRandom", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
        {"refreshOnMaxBasisSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
        {"oneBodyRdmRange", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
        {"printEigenvectorsRange",
         SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
        {"refreshIterations", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
        {"DoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"HPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"ParticleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"PPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
        {"SinglesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT());

IMPLEMENT_ALGORITHM(UnrestrictedEquationOfMotionSinglesFromRpa) {

  if (in.get<int64_t>("complexVersion", 1) == 1) {
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
  bool preconditionerRandom(in.get<int64_t>("preconditionerRandom", 0) == 1);
  double preconditionerRandomSigma(
      in.get<double>("preconditionerRandomSigma", 0.1));
  bool refreshOnMaxBasisSize(in.get<int64_t>("refreshOnMaxBasisSize", 0) == 1);
  std::vector<int> oneBodyRdmIndices(
      RangeParser(in.get<std::string>("oneBodyRdmRange", "")).getRange());
  int eigenStates(in.get<int64_t>("eigenstates", 1));
  const double energyConvergence(in.get<double>("energyConvergence", 1e-6)),
      amplitudesConvergence(in.get<double>("amplitudesConvergence", 1e-6));
  bool intermediates(in.get<int64_t>("intermediates", 1));
  unsigned int maxIterations(in.get<int64_t>("maxIterations", 32));
  unsigned int minIterations(in.get<int64_t>("minIterations", 1));
  std::vector<int> eigenvectorsIndices(
      RangeParser(in.get<std::string>("printEigenvectorsRange", ""))
          .getRange());
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
  LOG(0, "RPAEomDavid") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "RPAEomDavid") << "energyConvergence " << energyConvergence
                        << std::endl;
  LOG(0, "RPAEomDavid") << eigenStates << " eigen states" << std::endl;
  LOG(0, "RPAEomDavid") << "No: " << No << std::endl;
  LOG(0, "RPAEomDavid") << "Nv: " << Nv << std::endl;
  LOG(0, "RPAEomDavid") << "maxBasisSize: " << maxBasisSize << std::endl;

  // Get copy of couloumb integrals

  Tensor<double> *pVijab(in.get<Tensor<double> *>("HHPPCoulombIntegrals"));
  Tensor<F> cVijab(pVijab->order,
                   pVijab->lens,
                   pVijab->sym,
                   *Sisi4s::world,
                   pVijab->get_name());
  Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  // HF terms
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia(new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (isArgumentGiven("HPFockMatrix") && isArgumentGiven("HHFockMatrix")
      && isArgumentGiven("PPFockMatrix")) {
    LOG(0, "RPAEomDavid") << "Using non-canonical orbitals" << std::endl;

    Tensor<double> *realFia(in.get<Tensor<double> *>("HPFockMatrix"));
    Tensor<double> *realFab(in.get<Tensor<double> *>("PPFockMatrix"));
    Tensor<double> *realFij(in.get<Tensor<double> *>("HHFockMatrix"));
    toComplexTensor(*realFij, *Fij);
    toComplexTensor(*realFab, *Fab);
    toComplexTensor(*realFia, *Fia);
  } else {
    LOG(0, "RPAEomDavid") << "Using canonical orbitals" << std::endl;
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

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H.setFij(Fij)
      .setFab(Fab)
      .setFia(Fia)
      .setVijab(Vijab)
      .setTabij(&Tabij)
      .with_right_apply_intermediates(intermediates)
      .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::RPA);

  struct RpaH {
  public:
    SimilarityTransformedHamiltonian<F> *h;
    SFockVector<F> right_apply(SFockVector<F> &V) {
      return h->right_apply_hirata_RPA(V);
    }
  } rpaH;
  rpaH.h = &H;

  // INITIALIZE SIMILARITY PRECONDITIONER
  CcsdPreconditioner<F> P;
  P.setTai(&Tai)
      .setTabij(&Tabij)
      .setFij(Fij)
      .setFab(Fab)
      // Set coulomb integrals
      .setVijab(Vijab)
      // Set random information
      .setRandom(preconditionerRandom)
      .setRandomSigma(preconditionerRandomSigma);

  EigenSystemDavidsonMono<RpaH, CcsdPreconditioner<F>, SFockVector<F>>
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
    LOG(0, "RPAEomDavid") << "Refreshing on max basis size reaching"
                          << std::endl;
  }
  eigenSystem.run();

  if (eigenvectorsIndices.size() > 0) {

    for (auto &index : eigenvectorsIndices) {
      LOG(1, "RPAEomDavid") << "Writing out eigenvector " << index << std::endl;
      auto eigenState(eigenSystem.getRightEigenVectors()[index - 1]);
      TensorIo::writeText<F>("Rai-" + std::to_string(index) + ".tensor",
                             *eigenState.get(0),
                             "ij",
                             "",
                             " ");
    }
  }

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  NEW_FILE("EomRPAEnergies.dat") << "";
  for (auto &ev : eigenValues) {
    eigenCounter++;
    LOG(0, "RPAEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("EomRPAEnergies.dat")
        << eigenCounter << " " << ev.real() << " " << ev.imag() << std::endl;
  }
}
