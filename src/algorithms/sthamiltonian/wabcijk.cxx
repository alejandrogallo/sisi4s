#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <math/Complex.hpp>
#ifdef DEBUG
#define ST_DEBUG(msg) \
  LOG(1, "debug:STHam:") << __LINE__ << ":" << \
  "\x1b[33m" << msg << "\x1b[0m" << std::endl;
#else
#define ST_DEBUG(msg)
#endif

using namespace cc4s;

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getABCIJK() {
  if (Wabcijk) return Wabcijk;
  LOG(1, getAbbreviation()) << "Building Wabcijk" << std::endl;
  const int syms[]   = {NS,NS,NS,NS,NS,NS};
  const int vvvooo[] = {Nv,Nv,Nv,No,No,No};

  Wabcijk = NEW(CTF::Tensor<F>,  6, vvvooo, syms, *Cc4s::world, "Wabcijk");

  //Triples
  (*Wabcijk)["edfijk"] = 0.0;

  if (Fia) {
    ST_DEBUG("Fia * T3 * T1")
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defojk"] * (*Fia)["oh"] * (*Tai)["hi"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defoki"] * (*Fia)["oh"] * (*Tai)["hj"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defoij"] * (*Fia)["oh"] * (*Tai)["hk"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["hdeijk"] * (*Fia)["oh"] * (*Tai)["fo"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["hdfijk"] * (*Fia)["oh"] * (*Tai)["eo"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["hfeijk"] * (*Fia)["oh"] * (*Tai)["do"];
    ST_DEBUG("Fia * T2 * T2")
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["depk"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["dfpk"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["fepk"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["depj"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["dfpj"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["fepj"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["fepi"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["dfpi"] * (*Fia)["pg"];
    (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["depi"] * (*Fia)["pg"];
  }

  ST_DEBUG("T2 * V")
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Viajk)["ofij"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Viajk)["oeij"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Viajk)["odij"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Viajk)["ofik"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Viajk)["oeik"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Viajk)["odik"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Viajk)["ofkj"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Viajk)["oekj"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Viajk)["odkj"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdjk"] * (*Vabic)["efig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gejk"] * (*Vabic)["dfig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfjk"] * (*Vabic)["edig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdki"] * (*Vabic)["efjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geki"] * (*Vabic)["dfjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfki"] * (*Vabic)["edjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Vabic)["efkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Vabic)["dfkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Vabic)["edkg"];

  //Residum equations
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Fij)["oi"] * (*Tabcijk)["defojk"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Fij)["oj"] * (*Tabcijk)["defoik"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Fij)["ok"] * (*Tabcijk)["defoji"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Fab)["fg"] * (*Tabcijk)["gdeijk"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Fab)["eg"] * (*Tabcijk)["gdfijk"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Fab)["dg"] * (*Tabcijk)["gfeijk"];

  ST_DEBUG("T3 * V")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Vijkl)["opij"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Vijkl)["opik"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Vijkl)["opkj"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdepjk"] * (*Viajb)["pfig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Viajb)["peig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gfepjk"] * (*Viajb)["pdig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdepki"] * (*Viajb)["pfjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdfpki"] * (*Viajb)["pejg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gfepki"] * (*Viajb)["pdjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Viajb)["pfkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdfpij"] * (*Viajb)["pekg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gfepij"] * (*Viajb)["pdkg"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Vabcd)["efgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Vabcd)["dfgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Vabcd)["edgh"];

  ST_DEBUG("T2 * T1 * V")
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijkl)["opij"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijkl)["opij"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijkl)["opij"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijkl)["opik"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijkl)["opik"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijkl)["opik"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijkl)["opkj"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijkl)["opkj"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijkl)["opkj"];

  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Viajb)["ofih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hj"] * (*Viajb)["oeih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Tai)["hj"] * (*Viajb)["odih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Viajb)["ofih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hk"] * (*Viajb)["oeih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hk"] * (*Viajb)["odih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Viajb)["ofjh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hi"] * (*Viajb)["oejh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["hi"] * (*Viajb)["odjh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Viajb)["ofkh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hi"] * (*Viajb)["oekh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hi"] * (*Viajb)["odkh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Viajb)["ofjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hk"] * (*Viajb)["oejh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hk"] * (*Viajb)["odjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Viajb)["ofkh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hj"] * (*Viajb)["oekh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hj"] * (*Viajb)["odkh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["ep"] * (*Viajb)["pfig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["dp"] * (*Viajb)["pfig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["fp"] * (*Viajb)["peig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gejk"] * (*Tai)["fp"] * (*Viajb)["pdig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["dp"] * (*Viajb)["peig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["ep"] * (*Viajb)["pdig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["ep"] * (*Viajb)["pfjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["dp"] * (*Viajb)["pfjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdki"] * (*Tai)["fp"] * (*Viajb)["pejg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geki"] * (*Tai)["fp"] * (*Viajb)["pdjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfki"] * (*Tai)["dp"] * (*Viajb)["pejg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["ep"] * (*Viajb)["pdjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Viajb)["pfkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Viajb)["pfkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Viajb)["pekg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Viajb)["pdkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Viajb)["pekg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Viajb)["pdkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Vabcd)["efgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Vabcd)["dfgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Vabcd)["edgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Vabcd)["efgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Vabcd)["dfgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Vabcd)["edgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Vabcd)["efgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Vabcd)["dfgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Vabcd)["edgh"];

  ST_DEBUG("T3 * T1 * V")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Aj"] * (*Vijka)["opiA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ak"] * (*Vijka)["opiA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Ai"] * (*Vijka)["opjA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ai"] * (*Vijka)["opkA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Ak"] * (*Vijka)["opjA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Aj"] * (*Vijka)["opkA"];

  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tai)["fI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tai)["eI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tai)["dI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tai)["fI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tai)["eI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tai)["dI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["fI"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tai)["eI"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tai)["dI"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defojk"] * (*Tai)["hI"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defoki"] * (*Tai)["hI"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defoij"] * (*Tai)["hI"] * (*Vijka)["oIkh"];

  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tai)["Ai"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tai)["Ai"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tai)["Ai"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tai)["Aj"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tai)["Aj"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tai)["Aj"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["Ak"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tai)["Ak"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tai)["Ak"] * (*Viabc)["pdgA"];

  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["eI"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["dI"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["fI"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["fI"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["dI"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["eI"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdeijk"] * (*Tai)["hI"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdfijk"] * (*Tai)["hI"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gfeijk"] * (*Tai)["hI"] * (*Viabc)["Idgh"];

  ST_DEBUG("T2 * T2 * V")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gfjk"] * (*Tabij)["depI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gejk"] * (*Tabij)["dfpI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gdjk"] * (*Tabij)["fepI"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gfki"] * (*Tabij)["depI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["geki"] * (*Tabij)["dfpI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gdki"] * (*Tabij)["fepI"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Vijka)["pIkg"];

  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Vijka)["oIih"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Vijka)["oIjh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Vijka)["oIkh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Vijka)["oIkh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Vijka)["oIkh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Vijka)["oIkh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Vijka)["oIkh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Vijka)["oIkh"];

  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Viabc)["Idgh"];

  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Viabc)["ofhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfok"] * (*Tabij)["hAij"] * (*Viabc)["oehA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["feok"] * (*Tabij)["hAij"] * (*Viabc)["odhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Viabc)["ofhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["dfoj"] * (*Tabij)["hAik"] * (*Viabc)["oehA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["feoj"] * (*Tabij)["hAik"] * (*Viabc)["odhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Viabc)["ofhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"] * (*Viabc)["oehA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["feoi"] * (*Tabij)["hAkj"] * (*Viabc)["odhA"];
  ST_DEBUG("T3 * T2 * V (1)")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["Bfij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIk"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["Beij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIk"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["Bdij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIk"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["Bfik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIj"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["Beik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIj"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["Bdik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIj"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["Bfkj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIi"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["Bekj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIi"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["Bdkj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIi"];
  ST_DEBUG("T3 * T2 * V (2)")
  (*Wabcijk)["edfijk"] += ( + 0.25  ) * (*Tabcijk)["defopk"] * (*Tabij)["ABij"] * (*Vijab)["opAB"];
  (*Wabcijk)["edfijk"] += ( - 0.25  ) * (*Tabcijk)["defopj"] * (*Tabij)["ABik"] * (*Vijab)["opAB"];
  (*Wabcijk)["edfijk"] += ( - 0.25  ) * (*Tabcijk)["defopi"] * (*Tabij)["ABkj"] * (*Vijab)["opAB"];
  ST_DEBUG("T3 * T2 * V (3)")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["efJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIjk"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIjk"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["edJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIjk"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["efJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIki"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIki"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["edJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIki"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["efJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIij"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIij"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["edJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIij"];
  ST_DEBUG("T3 * T2 * V (4)")
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tabij)["AfJi"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpjk"] * (*Tabij)["AeJi"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepjk"] * (*Tabij)["AdJi"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepki"] * (*Tabij)["AfJj"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpki"] * (*Tabij)["AeJj"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepki"] * (*Tabij)["AdJj"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Tabij)["AfJk"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tabij)["AeJk"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["gfepij"] * (*Tabij)["AdJk"] * (*Vijab)["pJgA"];
  ST_DEBUG("T3 * T2 * V (5)")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defojk"] * (*Tabij)["hAJi"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defoki"] * (*Tabij)["hAJj"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defoij"] * (*Tabij)["hAJk"] * (*Vijab)["oJhA"];
  ST_DEBUG("T3 * T2 * V (6)")
  (*Wabcijk)["edfijk"] += ( + 0.25  ) * (*Tabcijk)["ghdijk"] * (*Tabij)["efIJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 0.25  ) * (*Tabcijk)["gheijk"] * (*Tabij)["dfIJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 0.25  ) * (*Tabcijk)["ghfijk"] * (*Tabij)["edIJ"] * (*Vijab)["IJgh"];
  ST_DEBUG("T3 * T2 * V (7)")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["gdeijk"] * (*Tabij)["hfIJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["gdfijk"] * (*Tabij)["heIJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["gfeijk"] * (*Tabij)["hdIJ"] * (*Vijab)["IJgh"];
  ST_DEBUG("T1 * T1 * T2 * Vhhhp")
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deIk"] * (*Tai)["gj"] * (*Tai)["fp"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfIk"] * (*Tai)["gj"] * (*Tai)["ep"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feIk"] * (*Tai)["gj"] * (*Tai)["dp"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deIj"] * (*Tai)["gk"] * (*Tai)["fp"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfIj"] * (*Tai)["gk"] * (*Tai)["ep"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feIj"] * (*Tai)["gk"] * (*Tai)["dp"] * (*Vijka)["pIig"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deIk"] * (*Tai)["gi"] * (*Tai)["fp"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfIk"] * (*Tai)["gi"] * (*Tai)["ep"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feIk"] * (*Tai)["gi"] * (*Tai)["dp"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deIj"] * (*Tai)["gi"] * (*Tai)["fp"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfIj"] * (*Tai)["gi"] * (*Tai)["ep"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feIj"] * (*Tai)["gi"] * (*Tai)["dp"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deIi"] * (*Tai)["gk"] * (*Tai)["fp"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfIi"] * (*Tai)["gk"] * (*Tai)["ep"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feIi"] * (*Tai)["gk"] * (*Tai)["dp"] * (*Vijka)["pIjg"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deIi"] * (*Tai)["gj"] * (*Tai)["fp"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfIi"] * (*Tai)["gj"] * (*Tai)["ep"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feIi"] * (*Tai)["gj"] * (*Tai)["dp"] * (*Vijka)["pIkg"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adjk"] * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["opiA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aejk"] * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["opiA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afjk"] * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["opiA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adki"] * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["opjA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aeki"] * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["opjA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afki"] * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["opjA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adij"] * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["opkA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aeij"] * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["opkA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afij"] * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["opkA"];
  ST_DEBUG("T1 * T1 * T2 * Vhppp")
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deIk"] * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfIk"] * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feIk"] * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deIj"] * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfIj"] * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feIj"] * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deIi"] * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Ifgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfIi"] * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Iegh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feIi"] * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Idgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Adjk"] * (*Tai)["gi"] * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Aejk"] * (*Tai)["gi"] * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adjk"] * (*Tai)["gi"] * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aejk"] * (*Tai)["gi"] * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Afjk"] * (*Tai)["gi"] * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afjk"] * (*Tai)["gi"] * (*Tai)["ep"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adik"] * (*Tai)["gj"] * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aeik"] * (*Tai)["gj"] * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Adik"] * (*Tai)["gj"] * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Aeik"] * (*Tai)["gj"] * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afik"] * (*Tai)["gj"] * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Afik"] * (*Tai)["gj"] * (*Tai)["ep"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Adji"] * (*Tai)["gk"] * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Aeji"] * (*Tai)["gk"] * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Adji"] * (*Tai)["gk"] * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Aeji"] * (*Tai)["gk"] * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Afji"] * (*Tai)["gk"] * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Afji"] * (*Tai)["gk"] * (*Tai)["ep"] * (*Viabc)["pdgA"];
  ST_DEBUG("T1 * T1 * T3 * V")
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defIJk"] * (*Tai)["gi"] * (*Tai)["hj"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["defIJj"] * (*Tai)["gi"] * (*Tai)["hk"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["defIJi"] * (*Tai)["gj"] * (*Tai)["hk"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["AdeJjk"] * (*Tai)["gi"] * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["AdfJjk"] * (*Tai)["gi"] * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["AfeJjk"] * (*Tai)["gi"] * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["AdeJik"] * (*Tai)["gj"] * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["AdfJik"] * (*Tai)["gj"] * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["AfeJik"] * (*Tai)["gj"] * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["AdeJji"] * (*Tai)["gk"] * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["AdfJji"] * (*Tai)["gk"] * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["AfeJji"] * (*Tai)["gk"] * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["defJjk"] * (*Tai)["gi"] * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defJik"] * (*Tai)["gj"] * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["defJji"] * (*Tai)["gk"] * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["ABdijk"] * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijab)["opAB"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabcijk)["ABeijk"] * (*Tai)["do"] * (*Tai)["fp"] * (*Vijab)["opAB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabcijk)["ABfijk"] * (*Tai)["do"] * (*Tai)["ep"] * (*Vijab)["opAB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabcijk)["Bdeijk"] * (*Tai)["fo"] * (*Tai)["hI"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["Bdfijk"] * (*Tai)["eo"] * (*Tai)["hI"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabcijk)["Bfeijk"] * (*Tai)["do"] * (*Tai)["hI"] * (*Vijab)["oIhB"];
  ST_DEBUG("T2 * T2 * T1 * V")
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gfjk"] * (*Tabij)["depI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gejk"] * (*Tabij)["dfpI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gdjk"] * (*Tabij)["fepI"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gfki"] * (*Tabij)["depI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["geki"] * (*Tabij)["dfpI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gdki"] * (*Tabij)["fepI"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Tai)["Bk"] * (*Vijab)["oIhB"];

  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];

  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["dfok"] * (*Tabij)["hAij"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["feok"] * (*Tabij)["hAij"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfoj"] * (*Tabij)["hAik"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["feoj"] * (*Tabij)["hAik"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 0.5  ) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"] * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( - 0.5  ) * (*Tabij)["feoi"] * (*Tabij)["hAkj"] * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["depk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["dfpk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["fepk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["depj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["dfpj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["fepj"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["fepi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["dfpi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["depi"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];

  ST_DEBUG("T1 * T1 * T1 * V")
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deJk"] * (*Tai)["hj"] * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfJk"] * (*Tai)["hj"] * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feJk"] * (*Tai)["hj"] * (*Tai)["dI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["deJj"] * (*Tai)["hk"] * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["dfJj"] * (*Tai)["hk"] * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["feJj"] * (*Tai)["hk"] * (*Tai)["dI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["deJi"] * (*Tai)["hk"] * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["dfJi"] * (*Tai)["hk"] * (*Tai)["eI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["feJi"] * (*Tai)["hk"] * (*Tai)["dI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Bdjk"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Bejk"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Bfjk"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Bdik"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Beik"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Bfik"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Bdji"] * (*Tai)["ep"] * (*Tai)["fI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( + 1.0  ) * (*Tabij)["Beji"] * (*Tai)["dp"] * (*Tai)["fI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];
  (*Wabcijk)["edfijk"] += ( - 1.0  ) * (*Tabij)["Bfji"] * (*Tai)["dp"] * (*Tai)["eI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];

  return Wabcijk;
}

// instantiate
template class SimilarityTransformedHamiltonian<cc4s::complex>;
template class SimilarityTransformedHamiltonian<double>;
