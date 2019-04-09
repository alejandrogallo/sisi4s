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
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdtEquationOfMotionDavidson);

CcsdtEquationOfMotionDavidson::CcsdtEquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
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
  bool preconditionerRandom(
    getIntegerArgument("preconditionerRandom", 0) == 1
  );
  double preconditionerRandomSigma(getRealArgument(
    "preconditionerRandomSigma", 0.1
  ));
  bool refreshOnMaxBasisSize(
    getIntegerArgument("refreshOnMaxBasisSize", 0) == 1
  );
  std::vector<int> oneBodyRdmIndices(
    RangeParser(getTextArgument("oneBodyRdmRange", "")).getRange()
  );
  int eigenStates(getIntegerArgument("eigenStates", 1));
  double ediff(getRealArgument("ediff", 1e-4));
  bool intermediates(getIntegerArgument("intermediates", 1));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));
  CTF::Tensor<double> *epsi(
    getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies")
  );
  CTF::Tensor<double> *epsa(
    getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies")
  );
  std::vector<int> refreshIterations(
    RangeParser(getTextArgument("refreshIterations", "")).getRange()
  );
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int  maxBasisSize(getIntegerArgument(
    "maxBasisSize", No*Nv + (No*(No - 1)/2 ) * (Nv * (Nv - 1)/2)
  ));

  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vv[] = {Nv, Nv};
  int ov[] = {No, Nv};
  int vo[] = {Nv,No};
  int oo[] = {No, No};
  int vvoo[] = {Nv,Nv,No,No};

  // Logging arguments
  LOG(0, "CcsdtEomDavid") << "Max iterations " << maxIterations << std::endl;
  LOG(0, "CcsdtEomDavid") << "ediff " << ediff << std::endl;
  LOG(0, "CcsdtEomDavid") << eigenStates << " eigen states" << std::endl;
  LOG(0, "CcsdtEomDavid") << "No: " << No << std::endl;
  LOG(0, "CcsdtEomDavid") << "Nv: " << Nv << std::endl;
  LOG(0, "CcsdtEomDavid") << "maxBasisSize: " << maxBasisSize << std::endl;


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


  // HF terms
  CTF::Tensor<F> *Fab(
    new CTF::Tensor<F>(2, vv, syms2, *Cc4s::world, "Fab")
  );
  CTF::Tensor<F> *Fij(
    new CTF::Tensor<F>(2, oo, syms2, *Cc4s::world, "Fij")
  );
  CTF::Tensor<F> *Fia(
    new CTF::Tensor<F>(2, ov, syms2, *Cc4s::world, "Fia")
  );

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    LOG(0, "CcsdtEomDavid") << "Using non-canonical orbitals" << std::endl;

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
    LOG(0, "CcsdtEomDavid") << "Using canonical orbitals" << std::endl;
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

  SimilarityTransformedHamiltonian<F> H(
    Fij, Fab, Fia,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci, NULL,
    intermediates
  );
  H.setTai(&Tai);
  H.setTabij(&Tabij);

  CcsdPreconditioner<F> P(
    Tai, Tabij, *Fij, *Fab, *Vabcd, *Viajb, *Vijab, *Vijkl
  );
  P.preconditionerRandom = preconditionerRandom;
  P.preconditionerRandomSigma = preconditionerRandomSigma;
  allocatedTensorArgument(
    "SinglesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(0))
  );
  allocatedTensorArgument(
    "DoublesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(1))
  );

  EigenSystemDavidsonMono<
    SimilarityTransformedHamiltonian<F>,
    CcsdPreconditioner<F>,
    CcsdtFockVector<F>
  > eigenSystem(
    &H,
    eigenStates,
    &P,
    ediff,
    maxBasisSize,
    maxIterations,
    minIterations
  );
  eigenSystem.refreshOnMaxBasisSize( refreshOnMaxBasisSize);
  if (eigenSystem.refreshOnMaxBasisSize()) {
    LOG(0, "CcsdtEomDavid") <<
      "Refreshing on max basis size reaching" << std::endl;
  }
  eigenSystem.run();

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  NEW_FILE("EomCcsdtEnergies.dat") << "";
  for (auto &ev: eigenValues) {
    eigenCounter++;
    LOG(0, "CcsdtEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
    FILE("EomCcsdtEnergies.dat") << eigenCounter <<
      " " << ev.real() << " " << ev.imag() << std::endl;
  }

}
