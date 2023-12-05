#include <algorithms/UCcsdIPEquationOfMotionDavidson.hpp>
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

ALGORITHM_REGISTRAR_DEFINITION(UCcsdIPEquationOfMotionDavidson);

IMPLEMENT_EMPTY_DRYRUN(UCcsdIPEquationOfMotionDavidson) {}

void UCcsdIPEquationOfMotionDavidson::run() {
  Data *Vabij(getArgumentData("HHPPCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  if (realVabij) {
    LOG(0, "IPEomDavid") << "Using real code" << std::endl;
    UCcsdIPEquationOfMotionDavidson::run<double>();
  } else {
    LOG(0, "IPEomDavid") << "Using complex code" << std::endl;
    UCcsdIPEquationOfMotionDavidson::run<complex>();
  }
}

template <typename F>
void UCcsdIPEquationOfMotionDavidson::run() {

  // Arguments
  bool refreshOnMaxBasisSize(getIntegerArgument("refreshOnMaxBasisSize", 0)
                             == 1);
  std::vector<int> oneBodyRdmIndices(
      RangeParser(getTextArgument("oneBodyRdmRange", "")).getRange());
  int eigenStates(getIntegerArgument("eigenstates", 1));
  bool intermediates(getIntegerArgument("intermediates", 1));
  const double energyConvergence(getRealArgument("energyConvergence", 1e-6)),
      amplitudesConvergence(getRealArgument("amplitudesConvergence", 1e-6));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));
  std::vector<int> eigenvectorsIndices(
      RangeParser(getTextArgument("printEigenvectorsRange", "")).getRange());
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
  LOG(0, "IPEomDavid") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "IPEomDavid") << "energyConvergence " << energyConvergence
                       << std::endl;
  LOG(0, "IPEomDavid") << eigenStates << " eigen states" << std::endl;
  LOG(0, "IPEomDavid") << "No: " << No << std::endl;
  LOG(0, "IPEomDavid") << "Nv: " << Nv << std::endl;
  LOG(0, "IPEomDavid") << "maxBasisSize: " << maxBasisSize << std::endl;

  // Get copy of couloumb integrals

  // Viabc
  Tensor<F> *Viabc = getTensorArgument<F, Tensor<F>>("HPPPCoulombIntegrals"),
            *Viajb = getTensorArgument<F, Tensor<F>>("HPHPCoulombIntegrals"),
            *Vaibc = getTensorArgument<F, Tensor<F>>("PHPPCoulombIntegrals"),
            *Viajk = getTensorArgument<F, Tensor<F>>("HPHHCoulombIntegrals"),
            *Vijab = getTensorArgument<F, Tensor<F>>("HHPPCoulombIntegrals"),
            *Vijka = getTensorArgument<F, Tensor<F>>("HHHPCoulombIntegrals"),
            *Vijkl = getTensorArgument<F, Tensor<F>>("HHHHCoulombIntegrals"),
            *Viabj = getTensorArgument<F, Tensor<F>>("HPPHCoulombIntegrals"),
            *Vaijb = getTensorArgument<F, Tensor<F>>("PHHPCoulombIntegrals"),
            // t
                *Tai = getTensorArgument<F, Tensor<F>>("SinglesAmplitudes"),
            *Tabij = getTensorArgument<F, Tensor<F>>("DoublesAmplitudes"),
            // HF terms
                *Fab = (new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab")),
            *Fij = (new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij")),
            *Fia = (new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (isArgumentGiven("HPFockMatrix") && isArgumentGiven("HHFockMatrix")
      && isArgumentGiven("PPFockMatrix")) {
    LOG(0, "IPEomDavid") << "Using non-canonical orbitals" << std::endl;

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
    LOG(0, "IPEomDavid") << "Using canonical orbitals" << std::endl;
    Fia = NULL;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

  H.setFij(Fij)
      .setFab(Fab)
      .setFia(Fia)
      .setViabc(Viabc)
      .setViabj(Viabj)
      .setViajk(Viajk)
      .setVijab(Vijab)
      .setVijka(Vijka)
      .setVijkl(Vijkl)
      // for intermediates
      .setViajb(Viajb)
      .setVaibc(Vaibc)
      .setVaijb(Vaijb)
      //
      .setTai(Tai)
      .setTabij(Tabij)
      // should we use intermediates of the Wabij etc?
      .with_right_apply_intermediates(intermediates)
      .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::CCSD);

  struct IPHamiltonian {
  public:
    SimilarityTransformedHamiltonian<F> *h;
    SDFockVector<F> right_apply(SDFockVector<F> &V) {
      return h->right_apply_CCSD_IP(V);
    }
  } ipH;
  ipH.h = &H;

  // INITIALIZE SIMILARITY PRECONDITIONER
  using _Preconditioner = IPCcsdPreconditioner<F>;
  _Preconditioner P;
  P.setTai(Tai)
      .setTabij(Tabij)
      .setFij(Fij)
      .setFab(Fab)
      // Set coulomb integrals
      .setVijab(Vijab);

  EigenSystemDavidsonMono<IPHamiltonian, _Preconditioner, SDFockVector<F>>
      eigenSystem(&ipH,
                  eigenStates,
                  &P,
                  amplitudesConvergence,
                  energyConvergence,
                  maxBasisSize,
                  maxIterations,
                  minIterations);
  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "IPEomDavid") << "Refreshing on max basis size reaching"
                         << std::endl;
  }
  eigenSystem.run();

  if (eigenvectorsIndices.size() > 0) {

    for (auto &index : eigenvectorsIndices) {
      LOG(1, "IPEomDavid") << "Writing out eigenvector " << index << std::endl;
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
  NEW_FILE("IPEomEnergies.dat") << "";
  for (auto &ev : eigenValues) {
    eigenCounter++;
    LOG(0, "IPEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("IPEomEnergies.dat")
        << eigenCounter << " " << ev.real() << " " << ev.imag() << std::endl;
  }
}
