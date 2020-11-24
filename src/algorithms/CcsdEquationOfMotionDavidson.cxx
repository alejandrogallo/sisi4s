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
  PTR(CTF::Tensor<F>) Vijkl, Vabcd, Vijka, Vijab, Viajk, Viajb, Viabc
                    , Vabic, Vabci, Vaibc, Vaibj, Viabj, Vijak, Vaijb
                    ;

  typedef struct { const std::string name; PTR(CTF::Tensor<F>) &data; } _Int;
  std::vector<_Int> requiredIntegrals =
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
                                          // structure factor
                                          , "CoulombVertex"
                                          , "GNorm"
                                          , "VG"
                                          , "structureFactorRange"
                                          , "structureFactorOnlySingles"
                                          , "structureFactorOnlyDoubles"
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

  const typename SimilarityTransformedHamiltonian<F>::StructureFactorSettings
    sfSettings = { getIntegerArgument("structureFactorOnlySingles", 0) == 1
                 , getIntegerArgument("structureFactorOnlyDoubles", 0) == 1
                 };

  const bool
    intermediates(getIntegerArgument("intermediates", 1))
  , refreshOnMaxBasisSize(getIntegerArgument("refreshOnMaxBasisSize", 0) == 1)
  , printEigenvectorsDoubles
      = getIntegerArgument("printEigenvectorsDoubles", 1) == 1
  ;

  const auto argToRange
    = [this](const std::string &a) -> std::vector<int> {
        return RangeParser(this->getTextArgument(a, "")).getRange();
      };

  const
  std::vector<int> refreshIterations(argToRange("refreshIterations"))
                 , oneBodyRdmIndices(argToRange("oneBodyRdmRange"))
                 , eigenvectorsIndices(argToRange("printEigenvectorsRange"))
                 , structureFactorIndices(argToRange("structureFactorRange"))
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
    , syms[] = {NS, NS, NS, NS}
    , vvoo[] = {Nv, Nv, No, No}
    , vv[] = {Nv, Nv}
    , ov[] = {No, Nv}
    , vo[] = {Nv, No}
    , oo[] = {No, No}
    ;

  const int maxBasisSize =
    getIntegerArgument("maxBasisSize", No*Nv
                                     + (No * (No - 1)/2)
                                     * (Nv * (Nv - 1)/2)
                                     );

  // Hilfsfunktion zum Rauschreiben von Tensoren
  const auto _writeText = TensorIo::writeText<F>;
  const auto _write_tensor =
    [&_writeText] (const std::string &name, char mode[], CTF::Tensor<F> &t) {
      _writeText(name, t, mode, "", " ");
    };



  EMIT() << YAML::Key << "No"            << YAML::Value << Nv
         << YAML::Key << "Nv"            << YAML::Value << No
         << YAML::Key << "maxBasisSize"  << YAML::Value << maxBasisSize
         << YAML::Key << "maxIterations" << YAML::Value << maxIterations
         << YAML::Key << "eigenStates"   << YAML::Value << eigenStates
         ;

  // Logging arguments
  LOGGER(0) << "max iter         : " << maxIterations << std::endl;
  LOGGER(0) << "energyConvergence: " << energyConvergence << std::endl;
  LOGGER(0) << "nroots           : " << eigenStates << std::endl;
  LOGGER(0) << "No               : " << No << std::endl;
  LOGGER(0) << "Nv               : " << Nv << std::endl;
  LOGGER(0) << "max basis        : " << maxBasisSize << std::endl;


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

  // INITIALIZE SIMILARITY TRANSFORMED HAMILTONIAN ===========================
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

  // INITIALIZE SIMILARITY PRECONDITIONER ====================================
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

  // INITIALIZE DAVIDSON SOLVER ==============================================
  EigenSystemDavidsonMono < SimilarityTransformedHamiltonian<F>
                          , CcsdPreconditioner<F>
                          , SDFockVector<F>
                          >
    eigenSystem( &H
               , eigenStates
               , &P
               , amplitudesConvergence
               , energyConvergence
               , maxBasisSize
               , maxIterations
               , minIterations
               )
               ;

  eigenSystem.refreshOnMaxBasisSize(refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize())
    LOGGER(0) << "Refreshing on max basis size reaching" << std::endl;

  eigenSystem.run();

  // 1-RDM ===================================================================
  if (oneBodyRdmIndices.size() > 0) {
    LOGGER(0) << "Calculating 1-RDM with left states "
              << " approximated by right" << std::endl;

    auto Ssquared(SzOperator<F>(No, Nv));

    _write_tensor("Szab.tensor",  "ij", *Ssquared.getAB());
    _write_tensor("Szij.tensor",  "ij", *Ssquared.getIJ());

    for (auto &index: oneBodyRdmIndices) {
      LOGGER(0) << "Calculating 1-RDM for state " << index << std::endl;

      const SDFockVector<F> *R(&eigenSystem.getRightEigenVectors()[index-1]);
      const SDFockVector<F> LApprox(R->conjugateTranspose());
      const SDFockVector<F> *L(&LApprox);

      EomOneBodyReducedDensityMatrix<F> Rho(&Tai, &Tabij, L, R);

      const std::string _name(std::to_string(index) + ".tensor");
      _write_tensor("Rhoia-" + _name, "ij", *Rho.getIA());
      _write_tensor("Rhoai-" + _name, "ij", *Rho.getAI());
      _write_tensor("Rhoab-" + _name, "ij", *Rho.getAB());
      _write_tensor("Rhoij-" + _name, "ij", *Rho.getAB());

      CTF::Scalar<F> s2;

      s2[""]  = (*Ssquared.getAB())["ab"] * (*Rho.getAB())["ba"];
      s2[""] += (*Ssquared.getIJ())["ij"] * (*Rho.getIJ())["ji"];

      F s2Val(s2.get_val());

      LOGGER(0) << "S^2 " << s2Val << std::endl;

    }
  }

  // WRITE AMPLITUDES ========================================================
  if (eigenvectorsIndices.size() > 0) {

    if (!printEigenvectorsDoubles) {
      LOGGER(0) << "Not writing out Rabij" << std::endl;
    }

    for (auto &index: eigenvectorsIndices) {
      LOGGER(1) << "Writing out eigenvector " << index << std::endl;
      auto eigenState(eigenSystem.getRightEigenVectors()[index-1]);

      _write_tensor( "Rai-" + std::to_string(index) + ".tensor"
                   , "ij"
                   , *eigenState.get(0)
                   );

      if (printEigenvectorsDoubles)
        _write_tensor( "Rabij-" + std::to_string(index) + ".tensor"
                     , "ijkl"
                     , *eigenState.get(1)
                     );

    }
  }

  // WRITE EIGENVALUES =======================================================
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


  // STRUCTURE FACTOR  =======================================================
  if (structureFactorIndices.size()) {
    auto GammaGqr(getTensorArgument<complex>("CoulombVertex"));
    complex madelungC;
    int64_t indices[] = {0};
    GammaGqr->read(1, indices, &madelungC);
    const double madelung(std::real(madelungC));
    LOGGER(1) << "Madelung : complex : " << madelungC << std::endl;
    LOGGER(1) << "Madelung : real    : " << madelung  << std::endl;
    const auto G(getTensorArgument<double>("GNorm"));

    // this is sqrt{ 1/kernel } without constants
    const double vMadelung = 1 / madelung;
    CTF::Tensor<double> V(G);
    V["G"] = (1.0 / 4.0 / M_PI) * (*G)["G"] * (*G)["G"];
    V.write(1, indices, &vMadelung);

    if (isArgumentGiven("VG")) allocatedTensorArgument("VG", &V);

    // set gamma in hamiltonian
    H.setGammaGqr(GammaGqr);
    for (auto &index: structureFactorIndices) {
      auto r(eigenSystem.getRightEigenVectors()[index-1]);
      LOGGER(1) << "Getting S for " << index << std::endl;
      auto S(H.structureFactor(r, sfSettings));
      LOGGER(1) << "Got S for " << index << std::endl;
      CTF::Scalar<F> energy;
      energy[""] = S.S["G"];
      LOGGER(1) << "Calculate Σ S[G]" << std::endl;
      const double value(std::real(energy.get_val()));
      const double braket(std::real(r.dot(r)));
      LOGGER(1) << "SF::energy "
                << "index"
                << " "
                << "S.energy"
                << " "
                << "value"
                << " "
                << "S.energy + value"
                << " "
                << "(S.energy + value) / r.dot(r)"
                << std::endl
                ;
      LOGGER(1) << "SF::energy "
                << index
                << " "
                << S.energy
                << " "
                << value
                << " "
                << S.energy + value
                << " "
                << (S.energy + value) / braket
                << std::endl
                ;
      // convert to real structure factor
      S.S["G"] = V["G"] * S.S["G"];
      _write_tensor("S-" + std::to_string(index) + ".tensor", "i", S.S);
    }
  }

}
