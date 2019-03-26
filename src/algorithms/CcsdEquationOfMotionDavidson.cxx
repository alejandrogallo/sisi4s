#include <algorithms/CcsdEquationOfMotionDavidson.hpp>
#include <algorithms/CcsdSimilarityTransformedHamiltonian.hpp>

#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/EigenSystemDavidson.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/MpiCommunicator.hpp>
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
using namespace tcc;

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
    //(*Fab)["aa"] = (*epsa)["a"];
    //(*Fij)["ii"] = (*epsi)["i"];
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
  // We initialize the T amplitudes here so that it is not necessary
  // to do a ccsd calculation before to do the CISD calculation.
  CTF::Tensor<F> Tai(2, vo, syms2, *Cc4s::world, "Tai");
  CTF::Tensor<F> Tabij(4, vvoo, syms4, *Cc4s::world, "Tabij");
  if (getIntegerArgument("CISD", 0) == 1) {
    LOG(0, "CcsdEomDavid") << "Calculating CISD" << std::endl;
    Tai["ai"] = 0.0;
    Tabij["abij"] = 0.0;
  } else {
    // Get the Uccsd amplitudes from the input file
    toComplexTensor(
      (*getTensorArgument<double, CTF::Tensor<double> >("SinglesAmplitudes")),
      Tai
    );
    toComplexTensor(
      (*getTensorArgument<double, CTF::Tensor<double> >("DoublesAmplitudes")),
      Tabij
    );
    //Tai["ai"] =
    //(*getTensorArgument<F, CTF::Tensor<F> >("SinglesAmplitudes"))["ai"];
    //Tabij["abij"] =
    //(*getTensorArgument<F, CTF::Tensor<F> >("DoublesAmplitudes"))["abij"];
  }

  bool intermediates(getIntegerArgument("intermediates", 1));
  CcsdSimilarityTransformedHamiltonian<F> H(
    &Tai, &Tabij, Fij, Fab, Fia,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci, NULL,
    intermediates
  );

  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));

  CcsdPreConditioner<F> P(
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
    CcsdPreConditioner<F>,
    FockVector<F>
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

      const FockVector<F> *R(&eigenSystem.getRightEigenVectors()[index-1]);
      const FockVector<F> LApprox(R->conjugateTranspose());
      const FockVector<F> *L(&LApprox);

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

template <typename F>
CcsdPreConditioner<F>::CcsdPreConditioner(
  CTF::Tensor<F> &Tai,
  CTF::Tensor<F> &Tabij,
  CTF::Tensor<F> &Fij,
  CTF::Tensor<F> &Fab,
  CTF::Tensor<F> &Vabcd,
  CTF::Tensor<F> &Viajb,
  CTF::Tensor<F> &Vijab,
  CTF::Tensor<F> &Vijkl
): diagonalH(
    std::vector<PTR(CTF::Tensor<F>)>(
      {NEW(CTF::Tensor<F>, Tai), NEW(CTF::Tensor<F>, Tabij)}
    ),
    std::vector<std::string>({"ai", "abij"})
  )
  {
  // pointers to singles and doubles tensors of diagonal part
  auto Dai( diagonalH.get(0) );
  auto Dabij( diagonalH.get(1) );

  // TODO: Maybe inster the Tai part to the diagonal

  // calculate diagonal elements of H
  (*Dai)["bi"] =  ( - 1.0 ) * Fij["ii"];
  (*Dai)["bi"] += ( + 1.0 ) * Fab["bb"];
/*  (*Dai)["bi"] += ( - 1.0 ) * Viajb["ibib"];
  (*Dai)["bi"] += ( + 1.0 ) * Tabij["cbli"] * Vijab["licb"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cblm"] * Vijab["lmcb"];
*/
  (*Dabij)["cdij"] =  ( - 1.0 ) * Fij["ii"];
  (*Dabij)["cdij"] += ( - 1.0 ) * Fij["jj"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Fab["cc"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Fab["dd"];
/*
  (*Dabij)["cdij"] += ( + 0.5 ) * Vijkl["ijij"];
  (*Dabij)["ccij"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["ccii"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["cdii"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( + 0.5 ) * Vabcd["cdcd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["efij"] * Vijab["ijef"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["ccij"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["ccii"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdii"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["cdmn"] * Vijab["mncd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecno"] * Vijab["noec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecno"] * Vijab["noec"];
*/
}


/**
 * \brief Comparator that should filter out zero values of the diagonal
 * matrix.
 * Zero values are treated as infinite so that they get appended to the
 * end of the list.
 */
template <typename F>
class EomDiagonalValueComparator {
public:
  bool operator ()(
    const std::pair<int, F> &a,
    const std::pair<int, F> &b
  ) {
    F A(
      std::abs(a.second) < 1E-13 ?
        std::numeric_limits<F>::infinity() : a.second
    );
    F B(
      std::abs(b.second) < 1E-13 ?
        std::numeric_limits<F>::infinity() : b.second
    );
    double diff(computeDifference(A, B));
    // maintain magnitude finite!
    double magnitude(std::abs(a.second)+std::abs(b.second));
    if (diff > +1E-13*magnitude) return true;
    if (diff < -1E-13*magnitude) return false;
    return a.first < b.first;
  }

  double computeDifference(const F &a, const F &b) { return b - a; }

};

//template<>
//class EomDiagonalValueComparator<double>;
//template<>
//class EomDiagonalValueComparator<complex>;

template<>
double EomDiagonalValueComparator<complex>::computeDifference(
    const complex &a,
    const complex &b
  ) {
  double diff(b.imag() + b.real() - a.imag() - a.real());
  return diff;
}



template <typename F>
std::vector<FockVector<F>>
CcsdPreConditioner<F>::getInitialBasis(const int eigenVectorsCount) {
  LOG(0, "CcsdPreConditioner") << "Getting initial basis " << std::endl;
  DefaultRandomEngine randomEngine;
  std::normal_distribution<double> normalDistribution(
    0.0, preconditionerRandomSigma
  );
  if (preconditionerRandom) {
    LOG(0, "CcsdPreConditioner") << "Randomizing initial guess" << std::endl;
  }
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    EomDiagonalValueComparator<F>()
  );
  int localElementsSize( localElements.size() );

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10*eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (
    int i(0);
    i < std::min(localElementsSize, trialEigenVectorsCount);
    ++i
  ) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (unsigned int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    EomDiagonalValueComparator<F>()
  );
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;

  int currentEigenVectorCount(0);
  unsigned int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t,F>> elements;
    if (communicator.getRank() == 0) {
      if ( b >= lowestElements.size() ) {
        throw EXCEPTION("No more elements to create initial basis");
      }
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    if (preconditionerRandom) {
      auto Rai(*basisElement.get(0));
      auto Rabij(*basisElement.get(1));
      setRandomTensor(Rai, normalDistribution, randomEngine);
      setRandomTensor(Rabij, normalDistribution, randomEngine);
      (*basisElement.get(0))["ai"] += Rai["ai"];
      (*basisElement.get(1))["abij"] += Rabij["abij"];
    }
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["abii"]=0.0;
    (*basisElement.get(1))["aaij"]=0.0;
    (*basisElement.get(1))["aaii"]=0.0;

    double preDot2(std::abs(basisElement.dot(basisElement)));
    // Antisymmetrize the new basis element
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["abji"];
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["baij"];

    OUT() << "\tnormPreSymmetrize=" << preDot2 << std::endl;

    double preDot3(std::abs(basisElement.dot(basisElement)));
    OUT() << "\tnormAfterSymmetrize=" << preDot3 << std::endl;

    OUT() << "\tbasisSize=" << basis.size() << std::endl;

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basis[j].dot(basisElement);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if ( std::abs(basisElementNorm) < 1e-10 ) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement = 1.0 / std::sqrt(basisElement.dot(basisElement))*basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if ( std::abs(basisElementNorm - double(1)) > 1e-10 * double(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);

  }

  // Now make sure that you are giving an antisymmetric basis
  // and also that it is grammschmited afterwards
  //LOG(1,"CcsdPreConditioner") << "Antisymmetrize basis" << std::endl;
  //for (unsigned int j(0); j < basis.size(); ++j) {
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["abji"];
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["baij"];
  //}

  //LOG(1,"CcsdPreConditioner") <<
      //"Performing Gramm Schmidt in the initial basis"
  //<< std::endl;
  //for (unsigned int b(0); b < basis.size(); ++b) {
    //V newVector(basis[b]);
    //for (unsigned int j(0); j < b; ++j) {
      //newVector -= basis[j] * basis[j].dot(basis[b]);
    //}
    //// normalize
    //basis[b] = 1 / std::sqrt(newVector.dot(newVector)) * newVector;
  //}

  return basis;
}

template <typename F>
FockVector<F>
CcsdPreConditioner<F>::getCorrection(
  const complex lambda, FockVector<F> &residuum
) {
  FockVector<F> w(diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
    public:
      DiagonalCorrection(const double lambda_): lambda(lambda_) {
      }
      F operator ()(const F residuumElement, const F diagonalElement) {
        return std::abs(lambda - diagonalElement) < 1E-4 ?
          0.0 : residuumElement / (lambda - diagonalElement);
      }
    protected:
      double lambda;
  } diagonalCorrection(std::real(lambda));

  FockVector<F> correction(diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.getComponentsCount(); ++c) {
    const char *indices( correction.componentIndices[c].c_str() );
    (*correction.get(c)).contract(
      1.0,
      *residuum.get(c),indices,
      *diagonalH.get(c),indices,
      0.0,indices,
      CTF::Bivar_Function<F>(diagonalCorrection)
    );
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["abii"]=0.0;
  (*correction.get(1))["aaij"]=0.0;
  (*correction.get(1))["aaii"]=0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["abij"] -= (*correction.get(1))["abji"];
  (*correction.get(1))["abij"] -= (*correction.get(1))["baij"];
  (*correction.get(1))["abij"] = 0.25 * (*correction.get(1))["abij"];

  return correction;
}


//// instantiate class
//template
//class CcsdPreConditioner<double>;

