#include <algorithms/CcsdEquationOfMotionDavidson.hpp>
#include <algorithms/CcsdSimilarityTransformedHamiltonian.hpp>
#include <algorithms/CcsdPreconditioner.hpp>

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

ALGORITHM_REGISTRAR_DEFINITION(CcsdEquationOfMotionDavidson);

CcsdEquationOfMotionDavidson::CcsdEquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
CcsdEquationOfMotionDavidson::~CcsdEquationOfMotionDavidson() {}


void CcsdEquationOfMotionDavidson::run() {

  if (getIntegerArgument("complexVersion", 1) == 1) {
    LOG(0, "CcsdEomDavid") << "Using complex code" << std::endl;
    CcsdEquationOfMotionDavidson::run<complex>();
  } else {
    LOG(0, "CcsdEomDavid") << "Using real code" << std::endl;
    CcsdEquationOfMotionDavidson::run<double>();
  }

}

template <typename F>
void CcsdEquationOfMotionDavidson::run() {

  // Get copy of couloumb integrals
  CTF::Tensor<double> *pVijkl(
    getTensorArgument<double, CTF::Tensor<double> >("HHHHCoulombIntegrals")
  );
  CTF::Tensor<F> cVijkl(
    pVijkl->order, pVijkl->lens, pVijkl->sym, *Cc4s::world,
    pVijkl->get_name()
  );
  CTF::Tensor<F> *Vijkl(&cVijkl);
  toComplexTensor(*pVijkl, *Vijkl);

  CTF::Tensor<double> *pVabcd(
    getTensorArgument<double, CTF::Tensor<double> >("PPPPCoulombIntegrals")
  );
  CTF::Tensor<F> cVabcd(
    pVabcd->order, pVabcd->lens, pVabcd->sym, *Cc4s::world,
    pVabcd->get_name()
  );
  CTF::Tensor<F> *Vabcd(&cVabcd);
  toComplexTensor(*pVabcd, *Vabcd);

  CTF::Tensor<double> *pVijka(
    getTensorArgument<double, CTF::Tensor<double> >("HHHPCoulombIntegrals")
  );
  CTF::Tensor<F> cVijka(
    pVijka->order, pVijka->lens, pVijka->sym, *Cc4s::world,
    pVijka->get_name()
  );
  CTF::Tensor<F> *Vijka(&cVijka);
  toComplexTensor(*pVijka, *Vijka);

  CTF::Tensor<double> *pVijab(
    getTensorArgument<double, CTF::Tensor<double> >("HHPPCoulombIntegrals")
  );
  CTF::Tensor<F> cVijab(
    pVijab->order, pVijab->lens, pVijab->sym, *Cc4s::world,
    pVijab->get_name()
  );
  CTF::Tensor<F> *Vijab(&cVijab);
  toComplexTensor(*pVijab, *Vijab);

  CTF::Tensor<double> *pViajk(
    getTensorArgument<double, CTF::Tensor<double> >("HPHHCoulombIntegrals")
  );
  CTF::Tensor<F> cViajk(
    pViajk->order, pViajk->lens, pViajk->sym, *Cc4s::world,
    pViajk->get_name()
  );
  CTF::Tensor<F> *Viajk(&cViajk);
  toComplexTensor(*pViajk, *Viajk);

  CTF::Tensor<double> *pViajb(
    getTensorArgument<double, CTF::Tensor<double> >("HPHPCoulombIntegrals")
  );
  CTF::Tensor<F> cViajb(
    pViajb->order, pViajb->lens, pViajb->sym, *Cc4s::world,
    pViajb->get_name()
  );
  CTF::Tensor<F> *Viajb(&cViajb);
  toComplexTensor(*pViajb, *Viajb);

  CTF::Tensor<double> *pViabc(
    getTensorArgument<double, CTF::Tensor<double> >("HPPPCoulombIntegrals")
  );
  CTF::Tensor<F> cViabc(
    pViabc->order, pViabc->lens, pViabc->sym, *Cc4s::world,
    pViabc->get_name()
  );
  CTF::Tensor<F> *Viabc(&cViabc);
  toComplexTensor(*pViabc, *Viabc);

  CTF::Tensor<double> *pVabic(
    getTensorArgument<double, CTF::Tensor<double> >("PPHPCoulombIntegrals")
  );
  CTF::Tensor<F> cVabic(
    pVabic->order, pVabic->lens, pVabic->sym, *Cc4s::world,
    pVabic->get_name()
  );
  CTF::Tensor<F> *Vabic(&cVabic);
  toComplexTensor(*pVabic, *Vabic);

  CTF::Tensor<double> *pVabci(
    getTensorArgument<double, CTF::Tensor<double> >("PPPHCoulombIntegrals")
  );
  CTF::Tensor<F> cVabci(
    pVabci->order, pVabci->lens, pVabci->sym, *Cc4s::world,
    pVabci->get_name()
  );
  CTF::Tensor<F> *Vabci(&cVabci);
  toComplexTensor(*pVabci, *Vabci);

  CTF::Tensor<double> *pVaibc(
    getTensorArgument<double, CTF::Tensor<double> >("PHPPCoulombIntegrals")
  );
  CTF::Tensor<F> cVaibc(
    pVaibc->order, pVaibc->lens, pVaibc->sym, *Cc4s::world,
    pVaibc->get_name()
  );
  CTF::Tensor<F> *Vaibc(&cVaibc);
  toComplexTensor(*pVaibc, *Vaibc);

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

  CTF::Tensor<double> *pVijak(
    getTensorArgument<double, CTF::Tensor<double> >("HHPHCoulombIntegrals")
  );
  CTF::Tensor<F> cVijak(
    pVijak->order, pVijak->lens, pVijak->sym, *Cc4s::world,
    pVijak->get_name()
  );
  CTF::Tensor<F> *Vijak(&cVijak);
  toComplexTensor(*pVijak, *Vijak);

  CTF::Tensor<double> *pVaijb(
    getTensorArgument<double, CTF::Tensor<double> >("PHHPCoulombIntegrals")
  );
  CTF::Tensor<F> cVaijb(
    pVaijb->order, pVaijb->lens, pVaijb->sym, *Cc4s::world,
    pVaijb->get_name()
  );
  CTF::Tensor<F> *Vaijb(&cVaijb);
  toComplexTensor(*pVaijb, *Vaijb);

  //CTF::Tensor<> *Vabij(
      //getTensorArgument<double, CTF::Tensor<>>("PPHHCoulombIntegrals"));

  // Get orbital energies
  CTF::Tensor<double> *epsi(
      getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies"));
  CTF::Tensor<double> *epsa(
      getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies"));
  int Nv(epsa->lens[0]), No(epsi->lens[0]);

  // HF terms
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int oo[] = {No, No};
  int kineticSyms[] = {NS, NS};
  CTF::Tensor<F> *Fab(
    new CTF::Tensor<F>(2, vv, kineticSyms, *Cc4s::world, "Fab")
  );
  CTF::Tensor<F> *Fij(
    new CTF::Tensor<F>(2, oo, kineticSyms, *Cc4s::world, "Fij")
  );
  CTF::Tensor<F> *Fia(
    new CTF::Tensor<F>(2, ov, kineticSyms, *Cc4s::world, "Fia")
  );

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    LOG(0, "CcsdEomDavid") << "Using non-canonical orbitals" << std::endl;

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
    LOG(0, "CcsdEomDavid") << "Using canonical orbitals" << std::endl;
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

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vo[] = {Nv,No};
  int vvoo[] = {Nv,Nv,No,No};
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

  bool intermediates(getIntegerArgument("intermediates", 1));
  CcsdSimilarityTransformedHamiltonian<F> H(
    Fij, Fab, Fia,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci, NULL,
    intermediates
  );
  H.setTai(&Tai);
  H.setTabij(&Tabij);

  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));

  CcsdPreconditioner<F> P(
    Tai, Tabij, *Fij, *Fab, *Vabcd, *Viajb, *Vijab, *Vijkl
  );
  P.preconditionerRandom = getIntegerArgument("preconditionerRandom", 0) == 1;
  P.preconditionerRandomSigma = getRealArgument(
    "preconditionerRandomSigma", 0.1
  );

  allocatedTensorArgument(
    "SinglesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(0))
  );
  allocatedTensorArgument(
    "DoublesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(1))
  );

  // Davidson solver
  int eigenStates(getIntegerArgument("eigenstates", 1));
  LOG(0, "CcsdEomDavid") << "Max iterations " << maxIterations << std::endl;
  double ediff(getRealArgument("ediff", 1e-4));
  LOG(0, "CcsdEomDavid") << "ediff " << ediff << std::endl;
  LOG(0, "CcsdEomDavid") << "Computing " << eigenStates << " eigen states"
                              << std::endl;
  std::vector<int> refreshIterations(
    RangeParser(getTextArgument("refreshIterations", "")).getRange()
  );
  EigenSystemDavidsonMono<
    CcsdSimilarityTransformedHamiltonian<F>,
    CcsdPreconditioner<F>,
    CcsdFockVector<F>
  > eigenSystem(
    &H,
    eigenStates,
    &P,
    ediff,
    getIntegerArgument(
      "maxBasisSize",
      No*Nv + (No*(No - 1)/2 ) * (Nv * (Nv - 1)/2)
    ),
    maxIterations,
    minIterations
  );
  eigenSystem.refreshOnMaxBasisSize(
      getIntegerArgument("refreshOnMaxBasisSize", 0) == 1
  );
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "CcsdEomDavid") <<
      "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();

  std::vector<int> oneBodyRdmIndices(
    RangeParser(getTextArgument("oneBodyRdmRange", "")).getRange()
  );

  if (oneBodyRdmIndices.size() > 0) {
    LOG(0, "CcsdEomDavid") << "Calculating 1-RDM with left states "
                           << " approximated by right" << std::endl;
    for (auto &index: oneBodyRdmIndices) {
      LOG(0, "CcsdEomDavid") << "Calculating 1-RDM for state " << index << std::endl;

      int syms[] = {NS, NS};
      CTF::Tensor<F> Rhoia(2, ov, syms, *Cc4s::world, "Rhoia");
      CTF::Tensor<F> Rhoai(2, vo, syms, *Cc4s::world, "Rhoai");
      CTF::Tensor<F> Rhoij(2, oo, syms, *Cc4s::world, "Rhoij");
      CTF::Tensor<F> Rhoab(2, vv, syms, *Cc4s::world, "Rhoab");

      const CcsdFockVector<F> *R(&eigenSystem.getRightEigenVectors()[index-1]);
      const CcsdFockVector<F> LApprox(R->conjugateTranspose());
      const CcsdFockVector<F> *L(&LApprox);

      Rhoia["ia"]  = 0;
      // this is 0 because r0 is 0
      //Rhoia["ia"] += (*L->get(0))["ia"];
      Rhoia["ia"] += (*L->get(1))["oifa"] * (*R->get(0))["fo"];
      // this is 0 because r0 is 0
      //Rhoia["ia"] += (*L->get(1))["oifa"] * Tai["fo"];
      //Rhoia.print(stdout);

      TensorIo::writeText<F>(
        "Rhoia-" + std::to_string(index) + ".tensor", Rhoia, "ij", "", " "
      );

      Rhoai["ai"]  = 0;
      Rhoai["ai"] += (*L->get(0))["ke"] * (*R->get(1))["eaki"];

      Rhoai["ai"] += (*L->get(0))["ke"] * (*R->get(0))["ak"] * Tai["ei"];
      Rhoai["ai"] += (*L->get(0))["ke"] * (*R->get(0))["ei"] * Tai["ak"];

      Rhoai["ai"] += (-0.5) * (*L->get(1))["kled"] * (*R->get(0))["di"] * Tabij["eakl"];
      Rhoai["ai"] += (-0.5) * (*L->get(1))["kled"] * (*R->get(0))["al"] * Tabij["edki"];

      Rhoai["ai"] += (-0.5) * (*L->get(1))["kled"] * Tai["di"] * (*R->get(1))["eakl"];
      Rhoai["ai"] += (-0.5) * (*L->get(1))["kled"] * Tai["al"] * (*R->get(1))["edki"];

      TensorIo::writeText<F>(
        "Rhoai-" + std::to_string(index) + ".tensor", Rhoai, "ij", "", " "
      );

      Rhoij["ij"]  = 0;
      Rhoij["ij"] += (*L->get(0))["je"] * (*R->get(0))["ei"];
      Rhoij["ij"] += 0.5 * (*L->get(1))["kjed"] * (*R->get(1))["edki"];
      Rhoij["ij"] += (*L->get(1))["kjed"] * (*R->get(0))["ek"] * Tai["di"];
      // This is not in the paper
      Rhoij["ij"] += (*L->get(1))["kjed"] * (*R->get(0))["di"] * Tai["ek"];

      TensorIo::writeText<F>(
        "Rhoij-" + std::to_string(index) + ".tensor", Rhoij, "ij", "", " "
      );

      Rhoab["ab"]  = 0;
      Rhoab["ab"] += (-1.0) * (*L->get(0))["ka"] * (*R->get(0))["bk"];
      Rhoab["ab"] += (-0.5) * (*L->get(1))["klea"] * (*R->get(1))["ebkl"];
      Rhoab["ab"] += (-1.0) * (*L->get(1))["klea"] * (*R->get(0))["ek"] * Tai["bl"];
      // This is not in the paper
      Rhoab["ab"] += (-1.0) * (*L->get(1))["klea"] * (*R->get(0))["bl"] * Tai["ek"];

      TensorIo::writeText<F>(
        "Rhoab-" + std::to_string(index) + ".tensor", Rhoab, "ij", "", " "
      );

    }
  }

  std::vector<int> eigenvectorsIndices(
    RangeParser(getTextArgument("printEigenvectorsRange", "")).getRange()
  );
  bool printEigenvectorsDoubles(getIntegerArgument("printEigenvectorsDoubles", 1) == 1);
  if (eigenvectorsIndices.size() > 0) {

    if (!printEigenvectorsDoubles) {
      LOG(0, "CcsdEomDavid") << "Not writing out Rabij" << std::endl;
    }

    for (auto &index: eigenvectorsIndices) {
      LOG(1, "CcsdEomDavid") << "Writing out eigenvector " << index << std::endl;
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
  NEW_FILE("EomCcsdEnergies.dat") << "";
  for (auto &ev: eigenValues) {
    eigenCounter++;
    LOG(0, "CcsdEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("EomCcsdEnergies.dat") << eigenCounter <<
      " " << ev.real() << " " << ev.imag() << std::endl;
  }

}
