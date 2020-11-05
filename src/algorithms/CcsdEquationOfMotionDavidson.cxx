#include <algorithms/CcsdEquationOfMotionDavidson.hpp>
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
#include <util/Emitter.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEquationOfMotionDavidson);
#define LOGGER(_l) LOG(_l, "CcsdEomDavid")

CcsdEquationOfMotionDavidson::CcsdEquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
CcsdEquationOfMotionDavidson::~CcsdEquationOfMotionDavidson() {}

template <typename F>
struct SpinOperator {
  SpinOperator(int No_, int Nv_): No(No_), Nv(Nv_) {};
  virtual PTR(CTF::Tensor<F>) getIJ() = 0;
  virtual PTR(CTF::Tensor<F>) getAB() = 0;
  PTR(CTF::Tensor<F>) Sab, Sij;
  int No, Nv;
};

template <typename F>
struct SzOperator: public SpinOperator<F> {
  SzOperator(int No, int Nv): SpinOperator<F>(No, Nv) {};
  PTR(CTF::Tensor<F>) getIJ() {
    if (this->Sij) return this->Sij;
    LOG(0, "SzOperator") << "Calculating Sz_ij" << std::endl;
    int oo[] = {this->No, this->No}, syms[] = {NS, NS};
    this->Sij = NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "Szij");
    (*this->Sij)["ii"] = 0.5;
    return this->Sij;
  }
  PTR(CTF::Tensor<F>) getAB() {
    if (this->Sab) return this->Sab;
    LOG(0, "SzOperator") << "Calculating Sz_ab" << std::endl;
    int vv[] = {this->Nv, this->Nv},  syms[] = {NS, NS};
    this->Sab = NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "Szab");
    (*this->Sab)["aa"] = 0.5;
    return this->Sab;
  }
};

void CcsdEquationOfMotionDavidson::run() {

  if (getIntegerArgument("complexVersion", 1) == 1) {
    LOGGER(0) << "Using complex code" << std::endl;
    CcsdEquationOfMotionDavidson::run<complex>();
  } else {
    LOGGER(0) << "Using real code" << std::endl;
    CcsdEquationOfMotionDavidson::run<double>();
  }

}

template <typename F>
void CcsdEquationOfMotionDavidson::run() {

  // initialize integrals, convert them to complex if we are using the
  // complex code
  PTR(CTF::Tensor<F>) Vijkl, Vabcd, Vijka, Vijab, Viajk, Viajb, Viabc,
                      Vabic, Vabci, Vaibc, Vaibj, Viabj, Vijak, Vaijb;

  typedef struct { const std::string name; PTR(CTF::Tensor<F>) &data; } _Int;
  std::vector<_Int>
  requiredIntegrals =
    { { "HHHHCoulombIntegrals", Vijkl }, { "PPPPCoulombIntegrals", Vabcd }
    , { "HHHPCoulombIntegrals", Vijka }, { "HHPPCoulombIntegrals", Vijab }
    , { "HPHHCoulombIntegrals", Viajk }, { "HPHPCoulombIntegrals", Viajb }
    , { "HPPPCoulombIntegrals", Viabc }, { "PPHPCoulombIntegrals", Vabic }
    , { "PPPHCoulombIntegrals", Vabci }, { "PHPPCoulombIntegrals", Vaibc }
    , { "PHPHCoulombIntegrals", Vaibj }, { "HPPHCoulombIntegrals", Viabj }
    , { "HHPHCoulombIntegrals", Vijak }, { "PHHPCoulombIntegrals", Vaijb }
    };

  std::vector<std::string> allArguments = { "complexVersion"
                                          , "oneBodyRdmRange"
                                          , "printEigenvectorsDoubles"
                                          , "printEigenvectorsRange"
                                          // Davidson solver
                                          , "amplitudesConvergence"
                                          , "energyConvergence"
                                          , "maxBasisSize"
                                          , "intermediates"
                                          , "eigenstates"
                                          , "refreshIterations"
                                          , "refreshOnMaxBasisSize"
                                          , "maxIterations"
                                          , "minIterations"
                                          // preconditioner
                                          , "preconditionerRandom"
                                          , "preconditionerRandomSigma"
                                          , "preconditionerSpinFlip"
                                          // T amplitudes
                                          , "SinglesAmplitudes"
                                          , "DoublesAmplitudes"
                                          // Fock Matrix
                                          , "ParticleEigenEnergies"
                                          , "HoleEigenEnergies"
                                          , "HPFockMatrix"
                                          , "HHFockMatrix"
                                          , "PPFockMatrix"
                                          };
  // possible integrals
  for (auto& i: requiredIntegrals) { allArguments.push_back(i.name); }
  checkArgumentsOrDie(allArguments);

  const struct { double sigma; bool random; bool spinFlip; }
  precSettings = { getRealArgument("preconditionerRandomSigma", 0.01)
                 , getIntegerArgument("preconditionerRandom", 0) == 1
                 , getIntegerArgument("preconditionerSpinFlip", 1) == 1
                 };

  const
  double energyConvergence(getRealArgument("energyConvergence", 1e-6))
       , amplitudesConvergence(getRealArgument("amplitudesConvergence", 1e-6))
       ;

  const bool
    intermediates(getIntegerArgument("intermediates", 1))
  , refreshOnMaxBasisSize(getIntegerArgument("refreshOnMaxBasisSize", 0) == 1)
  , printEigenvectorsDoubles(getIntegerArgument("printEigenvectorsDoubles", 1)
                             == 1)
  ;

  const std::function<std::vector<int>(const std::string)> argToRange
    = [&](const std::string a) {
        return RangeParser(getTextArgument(a, "")).getRange();
      };
  const
  std::vector<int> refreshIterations(argToRange("refreshIterations"))
                 , oneBodyRdmIndices(argToRange("oneBodyRdmRange"))
                 , eigenvectorsIndices(argToRange("printEigenvectorsRange"))
                 ;

  const unsigned
  int eigenStates(getIntegerArgument("eigenstates", 1))
    , maxIterations(getIntegerArgument("maxIterations", 32))
    , minIterations(getIntegerArgument("minIterations", 1))
    ;

  auto *epsi(getTensorArgument<double>("HoleEigenEnergies"))
     , *epsa(getTensorArgument<double>("ParticleEigenEnergies"))
     ;

  int Nv(epsa->lens[0])
    , No(epsi->lens[0])
    ;

  int syms[] = {NS, NS, NS, NS}
    , vvoo[] = {Nv, Nv, No, No}
    , vv[] = {Nv, Nv}
    , ov[] = {No, Nv}
    , vo[] = {Nv, No}
    , oo[] = {No, No}
    ;

  const int maxBasisSize =
    getIntegerArgument("maxBasisSize", No*Nv + (No * (No - 1)/2 ) *
                                               (Nv * (Nv - 1)/2)) ;

  EMIT() << YAML::Key << "No"            << YAML::Value << Nv
         << YAML::Key << "Nv"            << YAML::Value << No
         << YAML::Key << "maxBasisSize"  << YAML::Value << maxBasisSize
         << YAML::Key << "maxIterations" << YAML::Value << maxIterations
         << YAML::Key << "eigenStates"   << YAML::Value << eigenStates
         ;


  // Logging arguments
  LOGGER(0) << "max iter:  " << maxIterations << std::endl;
  LOGGER(0) << "energyConvergence:     " << energyConvergence << std::endl;
  LOGGER(0) << "nroots:    " << eigenStates << std::endl;
  LOGGER(0) << "No:        " << No << std::endl;
  LOGGER(0) << "Nv:        " << Nv << std::endl;
  LOGGER(0) << "max basis: " << maxBasisSize << std::endl;


  for (auto &integral: requiredIntegrals) {
    LOGGER(0) << "Converting " << integral.name << std::endl;
    auto in(getTensorArgument<double>(integral.name));
    integral.data = NEW(CTF::Tensor<F>, in->order, in->lens, in->sym,
                                        *in->wrld, in->get_name());
    toComplexTensor(*in, *integral.data);
  }

  // set up Fock matrix elements
  auto Fab(NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "Fab"));
  auto Fij(NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "Fij"));
  auto Fia(NEW(CTF::Tensor<F>, 2, ov, syms, *Cc4s::world, "Fia"));

  if (  isArgumentGiven("HPFockMatrix")
     && isArgumentGiven("HHFockMatrix")
     && isArgumentGiven("PPFockMatrix")
     ) {
    LOGGER(0) << "Using non-canonical orbitals" << std::endl;

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
    LOGGER(0) << "Using canonical orbitals" << std::endl;
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

  CTF::Tensor<F> Tai(2, vo, syms, *Cc4s::world, "Tai");
  CTF::Tensor<F> Tabij(4, vvoo, syms, *Cc4s::world, "Tabij");
  toComplexTensor(
    (*getTensorArgument<double, CTF::Tensor<double> >("SinglesAmplitudes")),
    Tai
  );
  toComplexTensor(
    (*getTensorArgument<double, CTF::Tensor<double> >("DoublesAmplitudes")),
    Tabij
  );

  // INITIALIZE SIMILARITY TRANSFORMED HAMILTONIAN
  SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);
  H
    // Set single particle integrals
    .setFij(Fij.get()).setFab(Fab.get()).setFia(Fia.get())
    // coulomb integrals setting
    .setVabcd(Vabcd.get()).setViajb(Viajb.get())
    .setVijab(Vijab.get()).setVijkl(Vijkl.get())
    .setVijka(Vijka.get()).setViabc(Viabc.get())
    .setViajk(Viajk.get()).setVabic(Vabic.get())
    .setVaibc(Vaibc.get()).setVaibj(Vaibj.get())
    .setViabj(Viabj.get()).setVijak(Vijak.get())
    .setVaijb(Vaijb.get()).setVabci(Vabci.get())
    // set dressing for the hamiltonian
    .setTai(&Tai).setTabij(&Tabij)
    // should we use intermediates of the Wabij etc?
    .setRightApplyIntermediates(intermediates)
    // Declare dressing of the hamiltonian so that we know that
    // Wai = Wabij = 0
    .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::CCSD)
  ;

  // INITIALIZE SIMILARITY PRECONDITIONER
  CcsdPreconditioner<F> P;
  P
    .setTai(&Tai).setTabij(&Tabij)
    .setFij(Fij.get()).setFab(Fab.get())

    // Set coulomb integrals
    //
    .setVabcd(Vabcd.get()).setViajb(Viajb.get())
    .setVijab(Vijab.get()).setVijkl(Vijkl.get())

    // preconditioner settings
    //
    .setRandom(precSettings.random)
    .setRandomSigma(precSettings.sigma)
    .setSpinFlip(precSettings.spinFlip)
  ;

  // INITIALIZE DAVIDSON SOLVER
  EigenSystemDavidsonMono < SimilarityTransformedHamiltonian<F>,
                            CcsdPreconditioner<F>,
                            SDFockVector<F> >
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
    LOGGER(0) << "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();


  if (oneBodyRdmIndices.size() > 0) {
    LOGGER(0) << "Calculating 1-RDM with left states "
              << " approximated by right" << std::endl;

    auto Ssquared(SzOperator<F>(No, Nv));
    TensorIo::writeText<F>("Szab.tensor", *Ssquared.getAB(), "ij", "", " ");
    TensorIo::writeText<F>("Szij.tensor", *Ssquared.getIJ(), "ij", "", " ");

    for (auto &index: oneBodyRdmIndices) {
      LOGGER(0) << "Calculating 1-RDM for state " << index << std::endl;

      const SDFockVector<F> *R(&eigenSystem.getRightEigenVectors()[index-1]);
      const SDFockVector<F> LApprox(R->conjugateTranspose());
      const SDFockVector<F> *L(&LApprox);

      EomOneBodyReducedDensityMatrix<F> Rho(&Tai, &Tabij, L, R);

      TensorIo::writeText<F>(
        "Rhoia-" + std::to_string(index) + ".tensor",
        *Rho.getIA(), "ij", "", " "
      );
      TensorIo::writeText<F>(
        "Rhoai-" + std::to_string(index) + ".tensor",
        *Rho.getAI(), "ij", "", " "
      );
      TensorIo::writeText<F>(
        "Rhoab-" + std::to_string(index) + ".tensor",
        *Rho.getAB(), "ij", "", " "
      );
      TensorIo::writeText<F>(
        "Rhoij-" + std::to_string(index) + ".tensor",
        *Rho.getAB(), "ij", "", " "
      );
      CTF::Scalar<F> s2;

      s2[""]  = (*Ssquared.getAB())["ab"] * (*Rho.getAB())["ba"];
      s2[""] += (*Ssquared.getIJ())["ij"] * (*Rho.getIJ())["ji"];
      F s2Val(s2.get_val());

      LOGGER(0) << "S^2 " << s2Val << std::endl;

    }
  }

  if (eigenvectorsIndices.size() > 0) {

    if (!printEigenvectorsDoubles) {
      LOGGER(0) << "Not writing out Rabij" << std::endl;
    }

    for (auto &index: eigenvectorsIndices) {
      LOGGER(1) << "Writing out eigenvector " << index << std::endl;
      auto eigenState(eigenSystem.getRightEigenVectors()[index-1]);
      TensorIo::writeText<F>(
        "Rai-" + std::to_string(index) + ".tensor",
        *eigenState.get(0),
        "ij", "", " "
      );
      if (printEigenvectorsDoubles) {
        TensorIo::writeText<F>(
          "Rabij-" + std::to_string(index) + ".tensor",
          *eigenState.get(1),
          "ijkl", "", " "
        );
      }
    }
  }

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  EMIT() << YAML::Key << "eigenValues"
         << YAML::Value << YAML::BeginSeq;
  for (auto &ev: eigenValues) {
    eigenCounter++;
    LOGGER(0) << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    EMIT() << YAML::BeginMap
              << YAML::Key << "real" << YAML::Value << ev.real()
              << YAML::Key << "imag" << YAML::Value << ev.imag()
           << YAML::EndMap;

  }
  EMIT() << YAML::EndSeq;

}
