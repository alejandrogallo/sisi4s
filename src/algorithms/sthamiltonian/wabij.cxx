#include <algorithms/SimilarityTransformedHamiltonian.hpp>
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
SimilarityTransformedHamiltonian<F>::getABIJ() {
  if (Wabij) return Wabij;

  LOG(1, getAbbreviation()) << "Building Wabij" << std::endl;

  Wabij = NEW(CTF::Tensor<F>, *Vabij);

  if (dressing == Dressing(CCSD)) {
    (*Wabij)["abij"] = 0.0;
    LOG(1, getAbbreviation()) << "Wabij = 0 since CCSD" << std::endl;
    return Wabij;
  }

  if (useStantonIntermediatesUCCSD()) {

    auto intermediates = getStantonIntermediatesUCCSD();
    (*Wabij)["abij"] = (*intermediates->getRabij())["abij"];
    //These are the residum equations
    (*Wabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
    (*Wabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
    (*Wabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
    (*Wabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];

  } else {

  // Equations are from hirata
  // WARNING: They are not optimized
  if (Fia) {
    (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Fia)["mf"] * (*Tai)["fi"];
    (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Fia)["mf"] * (*Tai)["fj"];
    (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["fcij"] * (*Fia)["mf"] * (*Tai)["dm"];
    (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["fdij"] * (*Fia)["mf"] * (*Tai)["cm"];
  }

  //These are the residum equations
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];

  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];

  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];

  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];

  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];

  (*Wabij)["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];

  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];

  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"] * (*Vijab)["opef"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdpj"] * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdpi"] * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"];

  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"] * (*Vijab)["mngh"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["hcij"] * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["hdij"] * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Tai)["dp"] * (*Vijab)["opef"];


  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Vabcd)["cdef"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];


  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Viajb)["ncje"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Viajb)["ncie"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["gcpi"] * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["npeg"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["gcpj"] * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["npeg"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["npeg"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["npeg"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];

  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Vabic)["cdie"];
  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Vabic)["cdje"];

  (*Wabij)["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajk)["mdij"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajk)["mcij"];

  // end without intermediates

  }

  if (Tabcijk) { // Begin Wabij with Triples Tabcijk

  if (Fia) {
    ST_DEBUG("Fia * T3")
    (*Wabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabcijk)["fcdmij"];
  }
  ST_DEBUG("T3 * Vijka")
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["ecdnoj"] * (*Vijka)["noie"];
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["ecdnoi"] * (*Vijka)["noje"];
  ST_DEBUG("T3 * Viabc")
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["efcoij"] * (*Viabc)["odef"];
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["efdoij"] * (*Viabc)["ocef"];

  ST_DEBUG("T3 * T1 * Vijab")
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["ecdnoj"] * (*Tai)["hi"] * (*Vijab)["noeh"];
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["ecdnoi"] * (*Tai)["hj"] * (*Vijab)["noeh"];
  (*Wabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["efcoij"] * (*Tai)["dp"] * (*Vijab)["opef"];
  (*Wabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["efdoij"] * (*Tai)["cp"] * (*Vijab)["opef"];
  (*Wabij)["cdij"] += ( + 1.0  ) * (*Tabcijk)["ecdnij"] * (*Tai)["gp"] * (*Vijab)["npeg"];

  } // end Wabij with Triples Tabcijk

  return Wabij;

}

// instantiate
template class SimilarityTransformedHamiltonian<complex>;
template class SimilarityTransformedHamiltonian<double>;
