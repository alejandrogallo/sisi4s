#include <equations/SimilarityTransformedHamiltonian.hpp>
#include <algorithms/CcsdPreconditioner.hpp>
#include <algorithms/OneBodyReducedDensityMatrix.hpp>
#include <algorithms/eom/CommonSpec.hpp>

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

namespace sisi4s {

enum EOM_TYPE { EA, EE, IP };

DEFSPEC(EOM, SPEC_IN(EOM_COMMON_SPEC_IN), SPEC_OUT());
DEFSTEP_METHODS_BEGIN(EOM)
template <typename F>
void run();
DEFSTEP_METHODS_END(EOM)

template <typename F>
void EOM::run() {

  // Arguments
  const std::string eom_type = in.get<std::string>("type");
  const bool restricted = in.get<bool>("restricted");
  bool refreshOnMaxBasisSize(in.get<bool>("refreshOnMaxBasisSize"));
  std::vector<int> oneBodyRdmIndices(
      RangeParser(in.get<std::string>("oneBodyRdmRange")).getRange());
  int eigenStates(in.get<int64_t>("eigenstates"));
  bool intermediates(in.get<bool>("intermediates"));
  const double energyConvergence(in.get<double>("energyConvergence")),
      amplitudesConvergence(in.get<double>("amplitudesConvergence"));
  unsigned int maxIterations(in.get<int64_t>("maxIterations"));
  unsigned int minIterations(in.get<int64_t>("minIterations"));
  std::vector<int> eigenvectorsIndices(
      RangeParser(in.get<std::string>("printEigenvectorsRange")).getRange());
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));
  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));
  std::vector<int> refreshIterations(
      RangeParser(in.get<std::string>("refreshIterations")).getRange());
  const int Nv(epsa->lens[0]), No(epsi->lens[0]);
  const int64_t maxBasisSize =
      in.get<int64_t>("maxBasisSize") == -1
          ? (No * Nv + (No * (No - 1) / 2) * (Nv * (Nv - 1) / 2))
          : in.get<int64_t>("maxBasisSize");

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int vo[] = {Nv, No};
  int oo[] = {No, No};
  int vvoo[] = {Nv, Nv, No, No};

  // Logging arguments
  LOG(0, "EOM") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "EOM") << "energyConvergence " << energyConvergence << std::endl;
  LOG(0, "EOM") << eigenStates << " eigen states" << std::endl;
  LOG(0, "EOM") << "No: " << No << std::endl;
  LOG(0, "EOM") << "Nv: " << Nv << std::endl;
  LOG(0, "EOM") << "maxBasisSize: " << maxBasisSize << std::endl;

  // Get copy of couloumb integrals

  // Viabc
  Tensor<F> *Viabc = in.get<Tensor<F> *>("HPPPCoulombIntegrals"),
            *Viajb = in.get<Tensor<F> *>("HPHPCoulombIntegrals"),
            *Vaibc = in.get<Tensor<F> *>("PHPPCoulombIntegrals"),
            *Vaibj = in.get<Tensor<F> *>("PHPHCoulombIntegrals"),
            *Viajk = in.get<Tensor<F> *>("HPHHCoulombIntegrals"),
            *Vijab = in.get<Tensor<F> *>("HHPPCoulombIntegrals"),
            *Vijka = in.get<Tensor<F> *>("HHHPCoulombIntegrals"),
            *Vijkl = in.get<Tensor<F> *>("HHHHCoulombIntegrals"),
            *Viabj = in.get<Tensor<F> *>("HPPHCoulombIntegrals"),
            *Vaijb = in.get<Tensor<F> *>("PHHPCoulombIntegrals"),
            *Vabci = in.get<Tensor<F> *>("PPPHCoulombIntegrals"),
            *Vabcd = in.get<Tensor<F> *>("PPPPCoulombIntegrals"),
            *Vijak = in.get<Tensor<F> *>("HHPHCoulombIntegrals"),
            // t
                *Tai = in.get<Tensor<F> *>("SinglesAmplitudes"),
            *Tabij = in.get<Tensor<F> *>("DoublesAmplitudes"),
            // HF terms
                *Fab = (new Tensor<F>(2, vv, syms2, *Sisi4s::world, "Fab")),
            *Fij = (new Tensor<F>(2, oo, syms2, *Sisi4s::world, "Fij")),
            *Fia = (new Tensor<F>(2, ov, syms2, *Sisi4s::world, "Fia"));

  if (in.present("HPFockMatrix") && in.present("HHFockMatrix")
      && in.present("PPFockMatrix")) {
    LOG(0, "EOM") << "Using non-canonical orbitals" << std::endl;

    Tensor<double> *realFia(in.get<Tensor<double> *>("HPFockMatrix"));
    Tensor<double> *realFab(in.get<Tensor<double> *>("PPFockMatrix"));
    Tensor<double> *realFij(in.get<Tensor<double> *>("HHFockMatrix"));
    toComplexTensor(*realFij, *Fij);
    toComplexTensor(*realFab, *Fab);
    toComplexTensor(*realFia, *Fia);
  } else {
    LOG(0, "EOM") << "Using canonical orbitals" << std::endl;
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
      .setVijak(Vijak)
      .setVijka(Vijka)
      .setVijkl(Vijkl)
      .setVaibj(Vaibj)
      .setVabcd(Vabcd)
      .setVabci(Vabci)
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

  if (in.present("TriplesAmplitudes")) {
    LOG(0, "EOM") << "Using CCSDT dressing" << std::endl;
    H
        // set CCSDT
        .setTabcijk(in.get<Tensor<F> *>("TriplesAmplitudes"))
        .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::CCSDT);
  }

#define MAKE_EOM_HAMILTONIAN(function_name, vector_type, var)                  \
  struct Eom_Hamiltonian {                                                     \
  public:                                                                      \
    SimilarityTransformedHamiltonian<F> *h;                                    \
    SDFockVector<F> right_apply(SDFockVector<F> &V) {                          \
      return h->function_name(V);                                              \
    }                                                                          \
  } var;                                                                       \
  var.h = &H;

#define RUN_EOM(eom_hamiltonian, vector_type)                                  \
  EigenSystemDavidsonMono<Eom_Hamiltonian, _Preconditioner, SDFockVector<F>>   \
      eigenSystem(&eom_hamiltonian,                                            \
                  eigenStates,                                                 \
                  &P,                                                          \
                  amplitudesConvergence,                                       \
                  energyConvergence,                                           \
                  maxBasisSize,                                                \
                  maxIterations,                                               \
                  minIterations);                                              \
  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);                    \
  if (eigenSystem.refreshOnMaxBasisSize()) {                                   \
    LOG(0, "EOM") << "Refreshing on max basis size reaching" << std::endl;     \
  }                                                                            \
  eigenSystem.run();

// INITIALIZE SIMILARITY PRECONDITIONER
#define MAKE_PRECONDITIONER(preconditioner_type)                               \
  using _Preconditioner = preconditioner_type;                                 \
  _Preconditioner P;                                                           \
  P.setTai(Tai).setTabij(Tabij).setFij(Fij).setFab(Fab).setVijab(Vijab);

  if (eom_type == "ip") {
    MAKE_EOM_HAMILTONIAN(right_apply_CCSD_IP, SDFockVector<F>, h);
    MAKE_PRECONDITIONER(IPCcsdPreconditioner<F>);
    RUN_EOM(h, SDFockVector<F>);
  } else if (eom_type == "ea") {
    MAKE_EOM_HAMILTONIAN(right_apply_CCSD_EA, SDFockVector<F>, h);
    MAKE_PRECONDITIONER(EACcsdPreconditioner<F>);
    RUN_EOM(h, SDFockVector<F>);
  } else if (eom_type == "ee") {
    MAKE_EOM_HAMILTONIAN(right_apply, SDFockVector<F>, h);
    MAKE_PRECONDITIONER(CcsdPreconditioner<F>);
    RUN_EOM(h, SDFockVector<F>);
  }
}

STEP_IMPLEMENT_RUN(EOM) {
  if (in.is_of_type<Tensor<double> *>("HHPPCoulombIntegrals")) {
    LOG(0, "EOM") << "Using real code" << std::endl;
    run<double>();
  } else {
    LOG(0, "EOM") << "Using complex code" << std::endl;
    run<complex>();
  }
}

} // namespace sisi4s
