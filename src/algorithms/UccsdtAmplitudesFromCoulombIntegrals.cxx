#include <algorithms/UccsdtAmplitudesFromCoulombIntegrals.hpp>
#include <unistd.h>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(UccsdtAmplitudesFromCoulombIntegrals);

UccsdtAmplitudesFromCoulombIntegrals::UccsdtAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesTriplesAlgorithm(argumentList) {
}

UccsdtAmplitudesFromCoulombIntegrals::~UccsdtAmplitudesFromCoulombIntegrals() {
}

void UccsdtAmplitudesFromCoulombIntegrals::run() {
  ClusterSinglesDoublesTriplesAlgorithm::run();
}

PTR(FockVector<complex>) UccsdtAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<complex>) &amplitudes
) {
  if (iterationStep == 0){
    LOG(1, getAbbreviation()) <<
       "WARNING: Using complex version of Uccsdt" << std::endl;
    LOG(1, getAbbreviation()) <<
       "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UccsdtAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<double>) &amplitudes
) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdtAmplitudesFromCoulombIntegrals::getResiduumTemplate(
  const int iterationStep, const PTR(const FockVector<F>) &amplitudes
) {
  // Equations from: hirata group
  // https://github.com/alejandrogallo/hirata

  CTF::Tensor<double> *epsi(
    getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies")
  );

  CTF::Tensor<double> *epsa(
    getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies")
  );

  // Get couloumb integrals
  auto Vijkl(getTensorArgument<F, CTF::Tensor<F> >("HHHHCoulombIntegrals"));
  auto Vabcd(getTensorArgument<F, CTF::Tensor<F> >("PPPPCoulombIntegrals"));
  auto Vijka(getTensorArgument<F, CTF::Tensor<F> >("HHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<F, CTF::Tensor<F> >("HHPPCoulombIntegrals"));
  auto Viajk(getTensorArgument<F, CTF::Tensor<F> >("HPHHCoulombIntegrals"));
  auto Viajb(getTensorArgument<F, CTF::Tensor<F> >("HPHPCoulombIntegrals"));
  auto Viabc(getTensorArgument<F, CTF::Tensor<F> >("HPPPCoulombIntegrals"));
  auto Vabij(getTensorArgument<F, CTF::Tensor<F> >("PPHHCoulombIntegrals"));
  auto Vabic(getTensorArgument<F, CTF::Tensor<F> >("PPHPCoulombIntegrals"));
  auto Viabj(getTensorArgument<F, CTF::Tensor<F> >("HPPHCoulombIntegrals"));
  auto Vaibc(getTensorArgument<F, CTF::Tensor<F> >("PHPPCoulombIntegrals"));
  auto Vijak(getTensorArgument<F, CTF::Tensor<F> >("HHPHCoulombIntegrals"));
  auto Vabci(getTensorArgument<F, CTF::Tensor<F> >("PPPHCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  CTF::Tensor<F> *Fab(
    new CTF::Tensor<F>(2, vv, syms, *Cc4s::world, "Fab")
  );
  CTF::Tensor<F> *Fij(
    new CTF::Tensor<F>(2, oo, syms, *Cc4s::world, "Fij")
  );
  CTF::Tensor<F> *Fia;

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    if (iterationStep == 0){
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    Fia = getTensorArgument<F, CTF::Tensor<F> >("HPFockMatrix");
    Fab = getTensorArgument<F, CTF::Tensor<F> >("PPFockMatrix");
    Fij = getTensorArgument<F, CTF::Tensor<F> >("HHFockMatrix");
  } else {
    if (iterationStep == 0){
      LOG(0, getAbbreviation()) << "Using hartree fock orbitals" << std::endl;
    }
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



  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  // amplitudes->get(0)
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");
  auto Tabcijk(amplitudes->get(2));
  Tabij->set_name("Tabcijk");


  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");
  auto Rabcijk(residuum->get(2));
  Rabcijk->set_name("Rabcijk");


  // Singles
  (*Rai)["bi"]  = 0.0;

  if (Fia) {
    (*Rai)["bi"] += ( + 1.0  ) * (*Fia)["ib"];
    (*Rai)["bi"] += ( + 1.0  ) * (*Fia)["kd"] * (*Tabij)["dbki"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
  }

  (*Rai)["bi"] += ( - 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];
  (*Rai)["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lbic"];
  (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
  (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
  (*Rai)["bi"] += ( + 0.25 ) * (*Tabcijk)["cdbmni"] * (*Vijab)["mncd"];
  (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
  (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
  (*Rai)["bi"] += ( - 0.5  ) * (*Tabij)["cblm"] * (*Tai)["fi"] * (*Vijab)["lmcf"];
  (*Rai)["bi"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tai)["bn"] * (*Vijab)["mncd"];
  (*Rai)["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
  (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Singles done" << std::endl;
#endif


  // Doubles
  (*Rabij)["cdij"]  = 0.0;

  if (Fia) {
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmj"] * (*Tai)["fi"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmi"] * (*Tai)["fj"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["fcij"] * (*Tai)["dm"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["fdij"] * (*Tai)["cm"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabcijk)["fcdmij"];
  }

  (*Rabij)["cdij"] += ( + 1.0  ) * (*Vabij)["cdij"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajk)["mdij"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajk)["mcij"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Vabic)["cdie"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Vabic)["cdje"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["ecdnoj"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["ecdnoi"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["efcoij"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["efdoij"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["ecdnoj"] * (*Tai)["hi"] * (*Vijab)["noeh"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["ecdnoi"] * (*Tai)["hj"] * (*Vijab)["noeh"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["efcoij"] * (*Tai)["dp"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["efdoij"] * (*Tai)["cp"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabcijk)["ecdnij"] * (*Tai)["gp"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Tabij)["gcpi"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tabij)["gdpi"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Tabij)["cdpj"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Tabij)["cdpi"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"] * (*Vijab)["mngh"];
  (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Tabij)["hcij"] * (*Vijab)["mofh"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Tabij)["hdij"] * (*Vijab)["mofh"];
  (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Tai)["dp"] * (*Vijab)["opef"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Doubles done" << std::endl;
#endif

  //Triples
  (*Rabcijk)["edfijk"] = 0.0;
  if (Fia) {
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defojk"] * (*Tai)["hi"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defoki"] * (*Tai)["hj"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defoij"] * (*Tai)["hk"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hdeijk"] * (*Tai)["fo"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hdfijk"] * (*Tai)["eo"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hfeijk"] * (*Tai)["do"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["depk"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["dfpk"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["fepk"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["depj"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["dfpj"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["fepj"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["fepi"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["dfpi"] * (*Fia)["pg"];
    (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["depi"] * (*Fia)["pg"];
  }
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deok"] * (*Viajk)["ofij"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfok"] * (*Viajk)["oeij"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feok"] * (*Viajk)["odij"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoj"] * (*Viajk)["ofik"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfoj"] * (*Viajk)["oeik"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feoj"] * (*Viajk)["odik"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoi"] * (*Viajk)["ofkj"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfoi"] * (*Viajk)["oekj"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feoi"] * (*Viajk)["odkj"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdjk"] * (*Vabic)["efig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gejk"] * (*Vabic)["dfig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfjk"] * (*Vabic)["edig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdki"] * (*Vabic)["efjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geki"] * (*Vabic)["dfjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfki"] * (*Vabic)["edjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdij"] * (*Vabic)["efkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Vabic)["dfkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfij"] * (*Vabic)["edkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fij)["oi"] * (*Tabcijk)["defojk"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Fij)["oj"] * (*Tabcijk)["defoik"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Fij)["ok"] * (*Tabcijk)["defoji"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Fab)["fg"] * (*Tabcijk)["gdeijk"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fab)["eg"] * (*Tabcijk)["gdfijk"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Fab)["dg"] * (*Tabcijk)["gfeijk"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Vijkl)["opij"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Vijkl)["opik"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Vijkl)["opkj"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdepjk"] * (*Viajb)["pfig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Viajb)["peig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gfepjk"] * (*Viajb)["pdig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdepki"] * (*Viajb)["pfjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdfpki"] * (*Viajb)["pejg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gfepki"] * (*Viajb)["pdjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Viajb)["pfkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdfpij"] * (*Viajb)["pekg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gfepij"] * (*Viajb)["pdkg"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Vabcd)["efgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Vabcd)["edgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijkl)["opij"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijkl)["opij"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijkl)["opij"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijkl)["opik"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijkl)["opik"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijkl)["opik"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijkl)["opkj"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijkl)["opkj"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijkl)["opkj"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Viajb)["ofih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hj"] * (*Viajb)["oeih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feok"] * (*Tai)["hj"] * (*Viajb)["odih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Viajb)["ofih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hk"] * (*Viajb)["oeih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hk"] * (*Viajb)["odih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Viajb)["ofjh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hi"] * (*Viajb)["oejh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["hi"] * (*Viajb)["odjh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Viajb)["ofkh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hi"] * (*Viajb)["oekh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hi"] * (*Viajb)["odkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Viajb)["ofjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hk"] * (*Viajb)["oejh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hk"] * (*Viajb)["odjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Viajb)["ofkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hj"] * (*Viajb)["oekh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hj"] * (*Viajb)["odkh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["ep"] * (*Viajb)["pfig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["dp"] * (*Viajb)["pfig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["fp"] * (*Viajb)["peig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gejk"] * (*Tai)["fp"] * (*Viajb)["pdig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["dp"] * (*Viajb)["peig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["ep"] * (*Viajb)["pdig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["ep"] * (*Viajb)["pfjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["dp"] * (*Viajb)["pfjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdki"] * (*Tai)["fp"] * (*Viajb)["pejg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geki"] * (*Tai)["fp"] * (*Viajb)["pdjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfki"] * (*Tai)["dp"] * (*Viajb)["pejg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["ep"] * (*Viajb)["pdjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Viajb)["pfkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Viajb)["pfkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Viajb)["pekg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Viajb)["pdkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Viajb)["pekg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Viajb)["pdkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Vabcd)["efgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Vabcd)["edgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Vabcd)["efgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Vabcd)["edgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Vabcd)["efgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Vabcd)["edgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Aj"] * (*Vijka)["opiA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ak"] * (*Vijka)["opiA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Ai"] * (*Vijka)["opjA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ai"] * (*Vijka)["opkA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Ak"] * (*Vijka)["opjA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Aj"] * (*Vijka)["opkA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tai)["fI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tai)["eI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tai)["dI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tai)["fI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tai)["eI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tai)["dI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["fI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tai)["eI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tai)["dI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["defojk"] * (*Tai)["hI"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["defoki"] * (*Tai)["hI"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["defoij"] * (*Tai)["hI"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tai)["Ai"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tai)["Ai"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tai)["Ai"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tai)["Aj"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tai)["Aj"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tai)["Aj"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["Ak"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tai)["Ak"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tai)["Ak"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["eI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["dI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["fI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["fI"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["dI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["eI"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdeijk"] * (*Tai)["hI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdfijk"] * (*Tai)["hI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gfeijk"] * (*Tai)["hI"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gfjk"] * (*Tabij)["depI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gejk"] * (*Tabij)["dfpI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gdjk"] * (*Tabij)["fepI"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gfki"] * (*Tabij)["depI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["geki"] * (*Tabij)["dfpI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gdki"] * (*Tabij)["fepI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Vijka)["oIkh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Viabc)["ofhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["dfok"] * (*Tabij)["hAij"] * (*Viabc)["oehA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["feok"] * (*Tabij)["hAij"] * (*Viabc)["odhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Viabc)["ofhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["dfoj"] * (*Tabij)["hAik"] * (*Viabc)["oehA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["feoj"] * (*Tabij)["hAik"] * (*Viabc)["odhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Viabc)["ofhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"] * (*Viabc)["oehA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["feoi"] * (*Tabij)["hAkj"] * (*Viabc)["odhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gdepIk"] * (*Tabij)["Bfij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gdfpIk"] * (*Tabij)["Beij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gfepIk"] * (*Tabij)["Bdij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gdepIj"] * (*Tabij)["Bfik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gdfpIj"] * (*Tabij)["Beik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gfepIj"] * (*Tabij)["Bdik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gdepIi"] * (*Tabij)["Bfkj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gdfpIi"] * (*Tabij)["Bekj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gfepIi"] * (*Tabij)["Bdkj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.25  ) * (*Tabcijk)["defopk"] * (*Tabij)["ABij"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( - 0.25  ) * (*Tabcijk)["defopj"] * (*Tabij)["ABik"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( - 0.25  ) * (*Tabcijk)["defopi"] * (*Tabij)["ABkj"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghdIjk"] * (*Tabij)["efJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gheIjk"] * (*Tabij)["dfJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghfIjk"] * (*Tabij)["edJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghdIki"] * (*Tabij)["efJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gheIki"] * (*Tabij)["dfJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghfIki"] * (*Tabij)["edJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["ghdIij"] * (*Tabij)["efJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gheIij"] * (*Tabij)["dfJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["ghfIij"] * (*Tabij)["edJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tabij)["AfJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tabij)["AeJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tabij)["AdJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tabij)["AfJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tabij)["AeJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tabij)["AdJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tabij)["AfJk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tabij)["AeJk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tabij)["AdJk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defojk"] * (*Tabij)["hAJi"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defoki"] * (*Tabij)["hAJj"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["defoij"] * (*Tabij)["hAJk"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.25  ) * (*Tabcijk)["ghdijk"] * (*Tabij)["efIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.25  ) * (*Tabcijk)["gheijk"] * (*Tabij)["dfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.25  ) * (*Tabcijk)["ghfijk"] * (*Tabij)["edIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabcijk)["gdeijk"] * (*Tabij)["hfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gdfijk"] * (*Tabij)["heIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabcijk)["gfeijk"] * (*Tabij)["hdIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["fp"] * (*Tabij)["deIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tabij)["dfIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tabij)["feIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["fp"] * (*Tabij)["deIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tabij)["dfIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tabij)["feIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["fp"] * (*Tabij)["deIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tabij)["dfIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tabij)["feIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["fp"] * (*Tabij)["deIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tabij)["dfIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tabij)["feIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["fp"] * (*Tabij)["deIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tabij)["dfIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tabij)["feIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["fp"] * (*Tabij)["deIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tabij)["dfIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tabij)["feIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adjk"] * (*Vijka)["opiA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aejk"] * (*Vijka)["opiA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afjk"] * (*Vijka)["opiA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adki"] * (*Vijka)["opjA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aeki"] * (*Vijka)["opjA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afki"] * (*Vijka)["opjA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adij"] * (*Vijka)["opkA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aeij"] * (*Vijka)["opkA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afij"] * (*Vijka)["opkA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["deIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["dfIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["feIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["deIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["dfIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["feIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["deIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["dfIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["feIi"] * (*Viabc)["Idgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tabij)["Adjk"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tabij)["Aejk"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["fp"] * (*Tabij)["Adjk"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["fp"] * (*Tabij)["Aejk"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tabij)["Afjk"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tabij)["Afjk"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tabij)["Adik"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tabij)["Aeik"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["fp"] * (*Tabij)["Adik"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["fp"] * (*Tabij)["Aeik"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tabij)["Afik"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tabij)["Afik"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tabij)["Adji"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tabij)["Aeji"] * (*Viabc)["pfgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["fp"] * (*Tabij)["Adji"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["fp"] * (*Tabij)["Aeji"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tabij)["Afji"] * (*Viabc)["pegA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tabij)["Afji"] * (*Viabc)["pdgA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabcijk)["defIJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabcijk)["defIJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabcijk)["defIJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["fp"] * (*Tabcijk)["AdeJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tabcijk)["AdfJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tabcijk)["AfeJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["fp"] * (*Tabcijk)["AdeJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tabcijk)["AdfJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tabcijk)["AfeJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["fp"] * (*Tabcijk)["AdeJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tabcijk)["AdfJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tabcijk)["AfeJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hI"] * (*Tabcijk)["defJjk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hI"] * (*Tabcijk)["defJik"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["hI"] * (*Tabcijk)["defJji"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabcijk)["ABdijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabcijk)["ABeijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabcijk)["ABfijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hI"] * (*Tabcijk)["Bdeijk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hI"] * (*Tabcijk)["Bdfijk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hI"] * (*Tabcijk)["Bfeijk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gfjk"] * (*Tabij)["depI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gejk"] * (*Tabij)["dfpI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gdjk"] * (*Tabij)["fepI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gfki"] * (*Tabij)["depI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["geki"] * (*Tabij)["dfpI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gdki"] * (*Tabij)["fepI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["dfok"] * (*Tabij)["hAij"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["feok"] * (*Tabij)["hAij"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["dfoj"] * (*Tabij)["hAik"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["feoj"] * (*Tabij)["hAik"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 0.5  ) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( - 0.5  ) * (*Tabij)["feoi"] * (*Tabij)["hAkj"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["depk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["dfpk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["fepk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["depj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["dfpj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["fepj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["fepi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["dfpi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["depi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["fI"] * (*Tabij)["deJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["eI"] * (*Tabij)["dfJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["dI"] * (*Tabij)["feJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["fI"] * (*Tabij)["deJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["eI"] * (*Tabij)["dfJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["dI"] * (*Tabij)["feJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["fI"] * (*Tabij)["deJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["eI"] * (*Tabij)["dfJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["dI"] * (*Tabij)["feJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tabij)["Bdjk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tabij)["Bejk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tabij)["Bfjk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tabij)["Bdik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tabij)["Beik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tabij)["Bfik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tabij)["Bdji"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( + 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tabij)["Beji"] * (*Vijab)["pIgB"];
  (*Rabcijk)["edfijk"] = ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tabij)["Bfji"] * (*Vijab)["pIgB"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Triples done" << std::endl;
#endif

  return residuum;
}



