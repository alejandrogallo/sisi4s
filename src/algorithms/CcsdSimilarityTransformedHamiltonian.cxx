/*
The equations in this file are taken from the following sources

[1] Isaiah Shavitt, Rodney J. Bartlett. Many-Body Methods in Chemistry and
    Physics: MBPT and Coupled-Cluster Theory. 2009
    PAGE: 439

[2] John F. Stanton, Rodney J. Bartlett. The equation of motion coupled‐cluster
    method. A systematic biorthogonal approach to molecular excitation
    energies, transition probabilities, and excited state properties. The
    Journal of Chemical Physics 7029--7039  1993
    TABLE 1
*/
#include <algorithms/CcsdSimilarityTransformedHamiltonian.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <tcc/DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>
#include <Cc4s.hpp>
#include <array>

#include <initializer_list>

using namespace CTF;
using namespace cc4s;


template <typename F>
CcsdSimilarityTransformedHamiltonian<F>::~CcsdSimilarityTransformedHamiltonian() {
}

template <typename F>
CcsdSimilarityTransformedHamiltonian<F>::CcsdSimilarityTransformedHamiltonian(
  CTF::Tensor<F> *Tai_,
  CTF::Tensor<F> *Tabij_,
  CTF::Tensor<F> *Fij_,
  CTF::Tensor<F> *Fab_,
  CTF::Tensor<F> *Fia_,
  CTF::Tensor<F> *Vabcd_,
  CTF::Tensor<F> *Viajb_,
  CTF::Tensor<F> *Vijab_,
  CTF::Tensor<F> *Vijkl_,
  CTF::Tensor<F> *Vijka_,
  CTF::Tensor<F> *Viabc_,
  CTF::Tensor<F> *Viajk_,
  CTF::Tensor<F> *Vabic_,
  CTF::Tensor<F> *Vaibc_,
  CTF::Tensor<F> *Vaibj_,
  CTF::Tensor<F> *Viabj_,
  CTF::Tensor<F> *Vijak_,
  CTF::Tensor<F> *Vaijb_,
  CTF::Tensor<F> *Vabci_,
  CTF::Tensor<F> *Vabij_,
  bool withIntermediates_,
  CcsdSimilarityTransformedHamiltonian::Dressing dressing_
):
  Tai(Tai_),
  Tabij(Tabij_),
  Fij(Fij_),
  Fab(Fab_),
  Fia(Fia_),
  Vabcd(Vabcd_),
  Viajb(Viajb_),
  Vijab(Vijab_),
  Vijkl(Vijkl_),
  Vijka(Vijka_),
  Viabc(Viabc_),
  Viajk(Viajk_),
  Vabic(Vabic_),
  Vaibc(Vaibc_),
  Vaibj(Vaibj_),
  Viabj(Viabj_),
  Vijak(Vijak_),
  Vaijb(Vaijb_),
  Vabci(Vabci_),
  Vabij(Vabij_),
  withIntermediates(withIntermediates_),
  dressing(dressing_)
{

  if (dressing == Dressing(CCSD)) {
    LOG(0, "CcsdSimilarityTransformedH") <<
      "Dressing is supposed to come from CCSD" << std::endl;
  } else if (dressing == Dressing(GENERAL)) {
    LOG(0, "CcsdSimilarityTransformedH") <<
      "Dressing is general Wai and Wabij are not 0" << std::endl;
  }

  No = Fij->lens[0];
  Nv = Fab->lens[0];

  Tau_abij = NEW(CTF::Tensor<F>, *Tabij);
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

}

template <typename F>
void CcsdSimilarityTransformedHamiltonian<F>::run() {
}

template <typename F>
CcsdFockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApply(
  CcsdFockVector<F> &R
) {
  return withIntermediates ? rightApplyIntermediates(R) : rightApplyHirata(R);
}

template <typename F>
CcsdFockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApplyHirata(
  CcsdFockVector<F> &R
) {
  CcsdFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  //checkAntisymmetry(*Rabij);

  // Contruct HR (one body part)
  // TODO: why "bi" not "ai"?
  (*HRai)["bi"]  = 0.0;

  // WIJ =====================================================================
  (*HRai)["bi"] += ( - 1.0 ) * (*Fij)["ki"] * (*Rai)["bk"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijka)["lmic"] * (*Rai)["bm"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Vijab)["mncd"] * (*Rai)["bn"];

  // WAB =====================================================================
  (*HRai)["bi"] += ( + 1.0 ) * (*Fab)["bc"] * (*Rai)["ci"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viabc)["lbce"] * (*Rai)["ei"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijab)["kmdf"] * (*Rai)["fi"];

  // WIABJ ===================================================================
  (*HRai)["bi"] += ( - 1.0 ) * (*Viajb)["kbid"] * (*Rai)["dk"];
  (*HRai)["bi"] += ( + 1.0 ) * (*Tabij)["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Vijka)["klie"] * (*Rai)["el"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Viabc)["lbce"] * (*Rai)["el"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Vijab)["lmcf"] * (*Rai)["fm"];

  // WIA =====================================================================
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijab)["lmcf"] * (*Rabij)["fbmi"];

  // WIJKA ===================================================================
  (*HRai)["bi"] += ( + 0.5 ) * (*Vijka)["klie"] * (*Rabij)["ebkl"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["ci"] * (*Vijab)["lmcf"] * (*Rabij)["fblm"];

  // WIABC ===================================================================
  (*HRai)["bi"] += ( + 0.5 ) * (*Viabc)["kbde"] * (*Rabij)["deki"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["bk"] * (*Vijab)["klef"] * (*Rabij)["efli"];


  //(*HRai)["ai"]  = 0.0; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  (*HRabij)["cdij"]  = 0.0;

  // Contruct HR (two body part)

  // WABCD ===================================================================
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vabcd)["cdef"] * (*Rabij)["efij"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["cm"] * (*Viabc)["mdfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["dm"] * (*Viabc)["mcfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["cdmn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];

  // WIJKL ===================================================================
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vijkl)["mnij"] * (*Rabij)["cdmn"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["efij"] * (*Vijab)["opef"] * (*Rabij)["cdop"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ej"] * (*Vijka)["noie"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["ei"] * (*Vijka)["noje"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["opef"] * (*Rabij)["cdop"];

  // WAB   ===================================================================
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fab)["de"] * (*Rabij)["ecij"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fab)["ce"] * (*Rabij)["edij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["ndeg"] * (*Rabij)["gcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["nceg"] * (*Rabij)["gdij"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edno"] * (*Vijab)["noeh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecno"] * (*Vijab)["noeh"] * (*Rabij)["hdij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hdij"];

  // WIJ   ===================================================================
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fij)["mi"] * (*Rabij)["cdmj"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fij)["mj"] * (*Rabij)["cdmi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["noie"] * (*Rabij)["cdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["noje"] * (*Rabij)["cdoi"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efoi"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efoj"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];

  // WIABJ ===================================================================
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mdif"] * (*Rabij)["fcmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mcif"] * (*Rabij)["fdmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mdjf"] * (*Rabij)["fcmi"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mcjf"] * (*Rabij)["fdmi"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnig"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnig"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnjg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnjg"] * (*Rabij)["gdni"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viabc)["ndeg"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viabc)["nceg"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viabc)["ndeg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viabc)["nceg"] * (*Rabij)["gdni"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];
  //--
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];

  //THREE_BODY_ONE ===========================================================
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rai)["ho"];

  //THREE_BODY_TWO ===========================================================
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edij"] * (*Vijab)["noeh"] * (*Rabij)["hcno"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecij"] * (*Vijab)["noeh"] * (*Rabij)["hdno"];

  //THREE_BODY_THREE =========================================================
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["cdmi"] * (*Vijka)["mnjg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fi"] * (*Vijab)["mofh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fj"] * (*Vijab)["mofh"] * (*Rai)["ho"];

  //THREE_BODY_FOUR ==========================================================
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mngh"] * (*Rabij)["ghnj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmj"] * (*Vijab)["mngh"] * (*Rabij)["ghni"];

  // WIAJK ===================================================================
  //--1
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajk)["mdij"] * (*Rai)["cm"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajk)["mcij"] * (*Rai)["dm"];
  //--2
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijkl)["mnij"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijkl)["mnij"] * (*Rai)["dn"];
  //--3
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viajb)["ndie"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viajb)["ncie"] * (*Rai)["dn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viajb)["ndje"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viajb)["ncje"] * (*Rai)["dn"];
  //--4
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijka)["noje"] * (*Rai)["do"];
  //--5
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["ocef"] * (*Rai)["do"];
  //--6
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
  //--7
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
  //--8
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["dp"];
  //--9
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  //--10
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
  //--11
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];

  // WABCI ===================================================================
  //--1
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Vabic)["cdie"] * (*Rai)["ej"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Vabic)["cdje"] * (*Rai)["ei"];
  //--2
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Vabcd)["cdef"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Vabcd)["cdef"] * (*Rai)["fi"];
  //--3
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdjf"] * (*Rai)["fi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcjf"] * (*Rai)["fi"];
  //--4
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  //--5
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  //--6
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  //--7
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  //--8
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hi"];
  //--9
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  //--10
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijab)["mngh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijab)["mngh"] * (*Rai)["hi"];
  //--11
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
   ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];

  // NON CANONICAL ORBITALS ==================================================

  if ( Fia ) {

    (*HRai)["bi"] +=  ( + 1.0  ) * (*Fia)["kd"] * (*Rabij)["dbki"];
    (*HRai)["bi"] +=  ( - 1.0  ) * (*Fia)["kd"] * (*Tai)["di"] * (*Rai)["bk"];
    (*HRai)["bi"] +=  ( - 1.0  ) * (*Fia)["kd"] * (*Tai)["bk"] * (*Rai)["di"];

    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tai)["fi"] * (*Rabij)["cdmj"];
    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tai)["fj"] * (*Rabij)["cdmi"];

    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tai)["dm"] * (*Rabij)["fcij"];
    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tai)["cm"] * (*Rabij)["fdij"];

    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["fdij"] * (*Rai)["cm"];
    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["fcij"] * (*Rai)["dm"];

    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmi"] * (*Rai)["fj"];
    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmj"] * (*Rai)["fi"];

  }

  return HR;

}

template <typename F>
CcsdFockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApplyIntermediates(
  CcsdFockVector<F> &R
) {
  CcsdFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  Wia = getIA();
  Wij = getIJ();
  Wab = getAB();
  Wabcd = getABCD();
  Wabci = getABCI();
  Waibc = getAIBC();
  Wiabj = getIABJ();
  Wiajk = getIAJK();
  Wijka = getIJKA();
  Wijkl = getIJKL();

  //checkAntisymmetry(*Rabij);

  (*HRai)["ai"]  = 0.0;
  (*HRai)["ai"] += (- 1.0) * (*Wij)["li"] * (*Rai)["al"];
  (*HRai)["ai"] += (*Wab)["ad"] * (*Rai)["di"];
  (*HRai)["ai"] += (*Wiabj)["ladi"] * (*Rai)["dl"];

  (*HRai)["ai"] += (*Wia)["ld"] * (*Rabij)["adil"];

  (*HRai)["ai"] += ( - 0.5 ) * (*Wijka)["lmid"] * (*Rabij)["adlm"];
  (*HRai)["ai"] += (   0.5 ) * (*Waibc)["alde"] * (*Rabij)["deil"];

  //(*HRai)["ai"]  = 0.0; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // 2 body part
  (*HRabij)["abij"]  = 0.0;

  // WABCD ===================================================================
  (*HRabij)["abij"] += (  0.5) * (*Wabcd)["abde"] * (*Rabij)["deij"];

  // WIJKL ===================================================================
  (*HRabij)["abij"] += (  0.5) * (*Wijkl)["lmij"] * (*Rabij)["ablm"];

  // WAB   ===================================================================
  (*HRabij)["abij"] += ( +1.0) * (*Wab)["bd"] * (*Rabij)["adij"];
  //P(ab)
  (*HRabij)["abij"] += ( -1.0) * (*Wab)["ad"] * (*Rabij)["bdij"];

  // WIJ   ===================================================================
  (*HRabij)["abij"] += ( -1.0) * (*Wij)["lj"] * (*Rabij)["abil"];
  //P(ij)
  (*HRabij)["abij"] +=           (*Wij)["li"] * (*Rabij)["abjl"];

  // WIABJ ===================================================================
  (*HRabij)["abij"] +=            (*Wiabj)["lbdj"] * (*Rabij)["adil"];
  //-P(ij)
  (*HRabij)["abij"] +=  ( -1.0) * (*Wiabj)["lbdi"] * (*Rabij)["adjl"];
  //-P(ab)
  (*HRabij)["abij"] +=  ( -1.0) * (*Wiabj)["ladj"] * (*Rabij)["bdil"];
  //P(ij)P(ab)
  (*HRabij)["abij"] +=            (*Wiabj)["ladi"] * (*Rabij)["bdjl"];

  //THREE_BODY_ONE ===========================================================
  (*HRabij)["abij"] +=
               (*Tabij)["afij"] * (*Rai)["em"] * (*Waibc)["bmfe"];
  // P(ab)
  (*HRabij)["abij"] +=
    ( -1.0) * (*Tabij)["bfij"] * (*Rai)["em"] * (*Waibc)["amfe"];

  //THREE_BODY_TWO ===========================================================
  (*HRabij)["abij"] +=
    ( -0.5) * (*Tabij)["fbij"] * (*Rabij)["eamn"] * (*Vijab)["nmfe"];
  // P(ab)
  (*HRabij)["abij"] +=
    ( +0.5) * (*Tabij)["faij"] * (*Rabij)["ebmn"] * (*Vijab)["nmfe"];

  //THREE_BODY_THREE =========================================================
  (*HRabij)["abij"] +=
    ( -1.0) * (*Tabij)["abin"] * (*Rai)["em"] * (*Wijka)["nmje"];
  // P(ij)
  (*HRabij)["abij"] +=
    ( +1.0) * (*Tabij)["abjn"] * (*Rai)["em"] * (*Wijka)["nmie"];

  //THREE_BODY_FOUR ==========================================================
  (*HRabij)["abij"] +=
    ( +0.5) * (*Tabij)["abjn"] * (*Vijab)["nmfe"] * (*Rabij)["feim"];
  // P(ij)
  (*HRabij)["abij"] +=
    ( -0.5) * (*Tabij)["abin"] * (*Vijab)["nmfe"] * (*Rabij)["fejm"];

  // WIAJK ===================================================================
  (*HRabij)["abij"] += (- 1.0 ) * (*Wiajk)["lbij"] * (*Rai)["al"];
  //P(ab)
  (*HRabij)["abij"] += (+ 1.0 ) * (*Wiajk)["laij"] * (*Rai)["bl"];

  // WABCI ===================================================================
  (*HRabij)["abij"] +=          (*Wabci)["abej"] * (*Rai)["ei"];
  //P(ij)
  (*HRabij)["abij"] += (-1.0) * (*Wabci)["abei"] * (*Rai)["ej"];

  return HR;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIJ() {
  if (Wij) return Wij;

  Wij = NEW(CTF::Tensor<F>, *Fij);

  if (true) {
    // This is the second row in diagram 10.55 in [1]
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wij from Wia" << std::endl;
    Wia = getIA();
    (*Wij)["ij"]  = (*Fij)["ij"];
    (*Wij)["ij"] += (*Vijka)["imje"] * (*Tai)["em"];
    (*Wij)["ij"] += (*Wia)["ie"] * (*Tai)["ej"];
    (*Wij)["ij"] += (  0.5) * (*Vijab)["imef"] * (*Tabij)["efjm"];
  } else {
    // This is the first row in diagram 10.55 in [1]
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wij" << std::endl;
    (*Wij)["ij"]  = (*Fij)["ij"];
    (*Wij)["ij"] += (*Vijka)["imje"] * (*Tai)["em"];
    if (Fia) {
      (*Wij)["ij"] += (*Fia)["ie"] * (*Tai)["ej"];
    }
    (*Wij)["ij"] += (  0.5) * (*Vijab)["imef"] * (*Tau_abij)["efjm"];
  }

  return Wij;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getAB() {
  if (Wab) return Wab;

  Wab   = NEW(CTF::Tensor<F>, *Fab);

  if (true) {
    //diagram (10.54) second line in [1]
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wab from Wia" << std::endl;
    Wia = getIA();
    (*Wab)["ab"]  = (*Fab)["ab"];
    (*Wab)["ab"] += (*Viabc)["mafb"] * (*Tai)["fm"];
    (*Wab)["ab"] += ( -1.0) * (*Wia)["mb"] * (*Tai)["am"];
    (*Wab)["ab"] += (- 0.5) * (*Vijab)["mnbe"] * (*Tabij)["aemn"];
  } else {
    //diagram (10.54) first line in [1]
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wab" << std::endl;
    (*Wab)["ab"]  = (*Fab)["ab"];
    (*Wab)["ab"] += (*Viabc)["mafb"] * (*Tai)["fm"];
    if (Fia) {
      (*Wab)["ab"] += ( -1.0) * (*Fia)["mb"] * (*Tai)["am"];
    }
    (*Wab)["ab"] += (- 0.5) * (*Vijab)["mnbe"] * (*Tau_abij)["aemn"];
  }

  return Wab;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getAI() {
  if (Wai) return Wai;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wai" << std::endl;

  int syms[] = {NS, NS};
  int ov[] = {Nv, No};
  CTF::Tensor<F> InitFai(2, ov, syms, *Cc4s::world, "InitFai");

  Wai = NEW(CTF::Tensor<F>, InitFai);

  (*Wai)["bi"] = 0.0;
  if (dressing == Dressing(CCSD)) {
    LOG(1, "CcsdSimilarityTransformedH") << "Wai = 0 since CCSD" << std::endl;
    return Wai;
  }

  // Equations from Hirata
  if (Fia) {
    (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["ib"];
    (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["kd"] * (*Tabij)["dbki"];
    (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
  }
  (*Wai)["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Wai)["bi"] += ( - 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];
  (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lbic"];
  (*Wai)["bi"] += ( + 0.5  ) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
  (*Wai)["bi"] += ( + 0.5  ) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
  (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
  (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
  (*Wai)["bi"] += ( - 0.5  ) * (*Tai)["fi"] * (*Tabij)["cblm"] * (*Vijab)["lmcf"];
  (*Wai)["bi"] += ( - 0.5  ) * (*Tai)["bn"] * (*Tabij)["cdmi"] * (*Vijab)["mncd"];
  (*Wai)["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
  (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];


  return Wai;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getABIJ() {
  if (Wabij) return Wabij;

  Wabij = NEW(CTF::Tensor<F>, *Vabij);

  if (dressing == Dressing(CCSD)) {
    (*Wabij)["abij"] = 0.0;
    LOG(1, "CcsdSimilarityTransformedH") <<
        "Wabij = 0 since CCSD" << std::endl;
    return Wabij;
  }

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

  return Wabij;

}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getABCIJK() {
  if (Wabcijk) return Wabcijk;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wabcijk from Wabci"
                                       << std::endl;
  int syms[] = {NS, NS, NS, NS};
  int vvvooo[] = {Nv,Nv,Nv,No,No,No};

  CTF::Tensor<F> Rabcijk(4, vvvooo, syms, *Cc4s::world, "Wabcijk");

  Wabcijk = NEW(CTF::Tensor<F>,  Rabcijk);
  Wabci = getABCI();

  (*Wabcijk)["abcijk"] = (+ 1.0) * (*Tabij)["aeij"] * (*Wabci)["bcek"];
  // we have to antisymmetrize a-b a-c and i-k j-k
  // a-b
  (*Wabcijk)["abcijk"] = (- 1.0) * (*Tabij)["beij"] * (*Wabci)["acek"];
  // a-c
  (*Wabcijk)["abcijk"] = (- 1.0) * (*Tabij)["ceij"] * (*Wabci)["baek"];

  return Wabcijk;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIA() {
  if (Wia) return Wia;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wia" << std::endl;

  int syms[] = {NS, NS};
  int ov[] = {No, Nv};
  CTF::Tensor<F> InitFia(2, ov, syms, *Cc4s::world, "InitFia");

  getABCIJK();

  Wia = NEW(CTF::Tensor<F>,  InitFia);

  //we need this one to construct the 2-body-amplitudes, not directly
  (*Wia)["ia"] = (*Vijab)["imae"] * (*Tai)["em"];
  if (Fia) {
    (*Wia)["ia"] += (*Fia)["ia"];
  }

  return Wia;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIJAB() {
  if (Wijab) return Wijab;

  LOG(0, "CcsdSimilarityTransformedH") << "Building Wijab = Vijab" << std::endl;
  Wijab = NEW(CTF::Tensor<F>, *Vijab);

  return Wijab;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getABCD() {
  if (Wabcd) return Wabcd;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wabcd" << std::endl;

  Wabcd = NEW(CTF::Tensor<F>, *Vabcd);
  // diagram 10.69 in [1]

  (*Wabcd)["abcd"]  = (*Vabcd)["abcd"];
  //-----------------------------------------------------------
  (*Wabcd)["abcd"] += (-1.0) * (*Vaibc)["amcd"] * (*Tai)["bm"];
  // P(ab)
  (*Wabcd)["abcd"] += ( 1.0) * (*Vaibc)["bmcd"] * (*Tai)["am"];
  //-----------------------------------------------------------
  (*Wabcd)["abcd"] += ( 0.5) * (*Vijab)["mncd"] * (*Tau_abij)["abmn"];

  return Wabcd;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getABCI() {
  if (Wabci) return Wabci;

  Wabci = NEW(CTF::Tensor<F>, *Vabci);

  bool wabciIntermediates(false);
  if (wabciIntermediates) {
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wabci from Wabcd and Wia"
                                         << std::endl;
    Wabcd = getABCD();
    Wia = getIA();
    //--1
    (*Wabci)["abci"]  = (*Vabci)["abci"];
    //--3
    (*Wabci)["abci"] += ( -1.0) * (*Vaibj)["amci"] * (*Tai)["bm"];
    (*Wabci)["abci"] += ( +1.0) * (*Vaibj)["bmci"] * (*Tai)["am"];
    //--6
    (*Wabci)["abci"] += ( +1.0) * (*Vaibc)["amce"] * (*Tabij)["ebmi"];
    (*Wabci)["abci"] += ( -1.0) * (*Vaibc)["bmce"] * (*Tabij)["eami"];
    //--9
    (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnce"] * (*Tai)["am"] * (*Tabij)["ebni"];
    (*Wabci)["abci"] += ( +1.0) * (*Vijab)["mnce"] * (*Tai)["bm"] * (*Tabij)["eani"];
    //--8
    (*Wabci)["abci"] += ( -1.0) * (*Wia)["mc"] * (*Tabij)["abmi"];
    //--2-4-10-11
    (*Wabci)["abci"] += ( +1.0) * (*Tai)["ei"] * (*Wabcd)["abce"];
    //--7-5
    (*Wabci)["abci"] += (  0.5 ) * (*Vijak)["nmci"] * (*Tau_abij)["abnm"];
  } else {
    LOG(0, "CcsdSimilarityTransformedH") << "Building Wabci" << std::endl;
    // from [1] first line of diagram 10.76, page 333
    //--1
    (*Wabci)["abci"]  = (*Vabci)["abci"];
    //--2
    (*Wabci)["abci"] += (*Vabcd)["abce"] * (*Tai)["ei"];
    //--3
    (*Wabci)["abci"] += ( -1.0) * (*Vaibj)["amci"] * (*Tai)["bm"];
    (*Wabci)["abci"] += ( +1.0) * (*Vaibj)["bmci"] * (*Tai)["am"];
    //--4
    (*Wabci)["abci"] += ( -1.0) * (*Vaibc)["amce"] * (*Tai)["bm"] * (*Tai)["ei"];
    (*Wabci)["abci"] += ( +1.0) * (*Vaibc)["bmce"] * (*Tai)["am"] * (*Tai)["ei"];
    //--5
    (*Wabci)["abci"] += ( +1.0) * (*Vijak)["mnci"] * (*Tai)["am"] * (*Tai)["bn"];
    //--5.1 (non canonical)
    if (Fia) {
      (*Wabci)["abci"] += ( -1.0) * (*Fia)["mc"] * (*Tabij)["abmi"];
    }
    //--6
    (*Wabci)["abci"] +=          (*Vaibc)["amce"] * (*Tabij)["ebmi"];
    (*Wabci)["abci"] += (-1.0) * (*Vaibc)["bmce"] * (*Tabij)["eami"];
    //--7
    (*Wabci)["abci"] += (  0.5) * (*Vijak)["mnci"] * (*Tabij)["abmn"];
    //--8
    (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnec"] * (*Tai)["em"] * (*Tabij)["abni"];
    //--9
    (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnce"] * (*Tai)["am"] * (*Tabij)["ebni"];
    (*Wabci)["abci"] += ( +1.0) * (*Vijab)["mnce"] * (*Tai)["bm"] * (*Tabij)["eani"];
    //--10
    (*Wabci)["abci"] += (  0.5) * (*Vijab)["mnce"] * (*Tai)["ei"] * (*Tabij)["abmn"];
    //--11
    (*Wabci)["abci"] +=           (*Vijab)["mnce"] * (*Tai)["am"] * (*Tai)["bn"] * (*Tai)["ei"];
  }

  return Wabci;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getAIBC() {
  if (Waibc) return Waibc;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Waibc" << std::endl;

  Waibc = NEW(CTF::Tensor<F>, *Vaibc);

  (*Waibc)["aibc"]  = (*Vaibc)["aibc"];
  (*Waibc)["aibc"] += ( -1.0) * (*Vijab)["mibc"] * (*Tai)["am"];

  return Waibc;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIABJ() {
  if (Wiabj) return Wiabj;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wiabj = Waijb" << std::endl;

  Wiabj = NEW(CTF::Tensor<F>, *Viabj);

  //[1] diagram (10.73)
  //This is not listed in the source book, however we can write it in terms
  //of Waijb since it should also have the simmetry of the Tabij amplitudes
  //and the Coulomb integrals Vpqrs
  //Taken directly from [2]
  (*Wiabj)["jabi"]  = (*Vaijb)["ajib"];
  (*Wiabj)["jabi"] += (*Vaibc)["ajeb"] * (*Tai)["ei"];
  (*Wiabj)["jabi"] += ( -1.0) * (*Vijka)["mjib"] * (*Tai)["am"];
  (*Wiabj)["jabi"] += ( -1.0) * (*Tai)["ei"] * (*Vijab)["mjeb"] * (*Tai)["am"];
  (*Wiabj)["jabi"] += ( -1.0) * (*Vijab)["mjeb"] * (*Tabij)["eaim"];

  return Wiabj;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIAJK() {
  if (Wiajk) return Wiajk;

  Wiajk = NEW(CTF::Tensor<F>, *Viajk);
  Wia = getIA();
  Wijkl = getIJKL();

  LOG(0, "CcsdSimilarityTransformedH") << "Building Wiajk from Wia and Wijkl" << std::endl;
  //This is built upon the already existing amplitudes
  //[1] diagram (10.79)
  //Takend directly from [2]
  //--1
  (*Wiajk)["iajk"]  = (*Viajk)["iajk"];
  //--6
  (*Wiajk)["iajk"] +=            (*Vijka)["imje"] * (*Tabij)["aekm"];
  (*Wiajk)["iajk"] += ( -1.0 ) * (*Vijka)["imke"] * (*Tabij)["aejm"];
  //    original
  //    (*Wiajk)["iajk"] +=            (*Vijka)["imje"] * (*Tabij)["aekm"];
  //    (*Wiajk)["iajk"] += ( -1.0 ) * (*Vijka)["jmie"] * (*Tabij)["aekm"];
  //--7-5
  (*Wiajk)["iajk"] += (  0.5 ) * (*Viabc)["iaef"] * (*Tau_abij)["efjk"];
  //--8
  (*Wiajk)["iajk"] += ( -1.0) * (*Wia)["ie"] * (*Tabij)["aejk"];
  //    original: (Problem, The diagram actually says that it
  //    should be Teajk and not Taejk, so that 'a' stays in the second
  //    vertex, so we have to either change a<>e or put a minus)
  //    (*Wiajk)["iajk"] += (*Wia)["ie"] * (*Tabij)["aejk"];
  //--2-4-10-11
  (*Wiajk)["iajk"] += (-1.0) * (*Tai)["am"] * (*Wijkl)["imjk"];
  //    original: (minus)
  //    (*Wiajk)["iajk"] += (*Tai)["am"] * (*Wijkl)["imjk"];
  //--3
  (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ek"] * (*Viajb)["iaje"];
  (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Viajb)["iake"];
  //     original:
  //     (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Viabj)["iaek"];
  //     (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ei"] * (*Viabj)["jaek"];
  //--9
  (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Tabij)["afmk"] * (*Vijab)["imef"];
  (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ek"] * (*Tabij)["afmj"] * (*Vijab)["imef"];
  //     original: Again it does not make any sense to do Pij, and the minus
  //     (*Wiajk)["iajk"] +=
  //       ( +1.0 ) * (*Tai)["ej"] * (*Tabij)["afmk"] * (*Vijab)["imef"];
  //     (*Wiajk)["iajk"] +=
  //       ( -1.0 ) * (*Tai)["ei"] * (*Tabij)["afmk"] * (*Vijab)["jmef"];

  return Wiajk;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIJKA() {
  if (Wijka) return Wijka;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wijka" << std::endl;

  Wijka = NEW(CTF::Tensor<F>, *Vijka);

  //Taken directly from[2]
  (*Wijka)["jkia"]  = (*Vijka)["jkia"];
  (*Wijka)["jkia"] += (*Tai)["ei"] * (*Vijab)["jkea"];

  return Wijka;
}

template <typename F>
PTR(CTF::Tensor<F>)
CcsdSimilarityTransformedHamiltonian<F>::getIJKL() {
  if (Wijkl) return Wijkl;
  LOG(0, "CcsdSimilarityTransformedH") << "Building Wijkl" << std::endl;

  Wijkl = NEW(CTF::Tensor<F>, *Vijkl);

  // diagram 10.70 in [1]
  //Taken directly from [2]
  (*Wijkl)["klij"]  = (*Vijkl)["klij"];
  //------------------------------------------------------------
  (*Wijkl)["klij"] +=           (*Tai)["ej"] * (*Vijka)["klie"];
  (*Wijkl)["klij"] += ( -1.0) * (*Tai)["ei"] * (*Vijka)["klje"];
  //------------------------------------------------------------
  (*Wijkl)["klij"] += ( 0.5 ) * (*Tau_abij)["efij"] * (*Vijab)["klef"];

  return Wijkl;
}

template <typename F>
CcsdFockVector<F> CcsdSimilarityTransformedHamiltonian<F>::leftApplyHirata(
  CcsdFockVector<F> &L
) {
  CcsdFockVector<F> LH(L);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Lia( L.get(0) );
  PTR(CTF::Tensor<F>) Lijab( L.get(1) );
  PTR(CTF::Tensor<F>) LHia( LH.get(0) );
  PTR(CTF::Tensor<F>) LHijab( LH.get(1) );

  // Contruct HR (one body part)
  (*LHia)["ja"]  = 0.0;
  (*LHia)["ja"] += ( - 1.0  ) * (*Fij)["jk"] * (*Lia)["ka"];
  (*LHia)["ja"] += ( + 1.0  ) * (*Fab)["ca"] * (*Lia)["jc"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Viajb)["jcla"] * (*Lia)["lc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Viajk)["jclm"] * (*Lijab)["mlca"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Vabic)["cdma"] * (*Lijab)["mjdc"];
  (*LHia)["ja"] += ( + 1.0  ) * (*Tabij)["cdmn"] * (*Vijab)["njda"] * (*Lia)["mc"];
  (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["njcd"] * (*Lia)["ma"];
  (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["mnda"] * (*Lia)["jc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["njoa"] * (*Lijab)["omdc"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Vijka)["njod"] * (*Lijab)["omca"];
  (*LHia)["ja"] += ( - 0.25 ) * (*Tabij)["cdmn"] * (*Vijka)["mnoa"] * (*Lijab)["ojdc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["mnod"] * (*Lijab)["ojca"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["jgda"] * (*Lijab)["nmgc"];
  (*LHia)["ja"] += ( - 0.25 ) * (*Tabij)["cdmn"] * (*Viabc)["jgcd"] * (*Lijab)["nmga"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Viabc)["ngda"] * (*Lijab)["mjgc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["ngcd"] * (*Lijab)["mjga"];

  // Contruct HR (two body part)
  (*LHijab)["klab"]  = 0.0;
  (*LHijab)["klab"] += ( - 1.0  ) * (*Vijka)["klmb"] * (*Lia)["ma"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Vijka)["klma"] * (*Lia)["mb"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viabc)["keab"] * (*Lia)["le"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viabc)["leab"] * (*Lia)["ke"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Fij)["km"] * (*Lijab)["mlab"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Fij)["lm"] * (*Lijab)["mkab"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Fab)["eb"] * (*Lijab)["klea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Fab)["ea"] * (*Lijab)["kleb"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Vijkl)["klmn"] * (*Lijab)["nmab"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viajb)["kenb"] * (*Lijab)["nlea"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viajb)["kena"] * (*Lijab)["nleb"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viajb)["lenb"] * (*Lijab)["nkea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viajb)["lena"] * (*Lijab)["nkeb"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Vabcd)["efab"] * (*Lijab)["klfe"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfb"] * (*Lijab)["poea"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfa"] * (*Lijab)["poeb"];
  (*LHijab)["klab"] += ( - 0.25 ) * (*Tabij)["efop"] * (*Vijab)["klef"] * (*Lijab)["poab"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkab"] * (*Lijab)["olfe"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plab"] * (*Lijab)["okfe"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfb"] * (*Lijab)["olea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfa"] * (*Lijab)["oleb"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfb"] * (*Lijab)["okea"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfa"] * (*Lijab)["okeb"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkef"] * (*Lijab)["olab"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plef"] * (*Lijab)["okab"];
  (*LHijab)["klab"] += ( - 0.25 ) * (*Tabij)["efop"] * (*Vijab)["opab"] * (*Lijab)["klfe"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfb"] * (*Lijab)["klea"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfa"] * (*Lijab)["kleb"];

  return LH;
}

template <typename F>
CcsdtFockVector<F> CcsdtSimilarityTransformedHamiltonian<F>::rightApply(
  CcsdtFockVector<F> &R
) {
  return rightApplyHirata(R);
}

template <typename F>
CcsdtFockVector<F> CcsdtSimilarityTransformedHamiltonian<F>::rightApplyHirata(
  CcsdtFockVector<F> &R
) {

  //: BEGIN SINGLES
  ( + 0.25  ) * Vijab["klef"] * Rabcijk["feblki"];

  //: BEGIN DOUBLES
  ( + 1.0  ) * Fia["mf"] * Rabcijk["fcdmij"];
  ( + 1.0  ) * Fab["fg"] * Rabcijk["gdeijk"];
  ( - 1.0  ) * Fab["eg"] * Rabcijk["gdfijk"];
  ( - 1.0  ) * Fab["dg"] * Rabcijk["gfeijk"];
  ( - 1.0  ) * Fij["oi"] * Rabcijk["defojk"];
  ( + 1.0  ) * Fij["oj"] * Rabcijk["defoik"];
  ( + 1.0  ) * Fij["ok"] * Rabcijk["defoji"];

  ( - 0.5  ) * Vijka["mnig"] * Rabcijk["gcdnmj"];
  ( + 0.5  ) * Vijka["mnjg"] * Rabcijk["gcdnmi"];

  ( + 0.5  ) * Viabc["mdfg"] * Rabcijk["gfcmij"];
  ( - 0.5  ) * Viabc["mcfg"] * Rabcijk["gfdmij"];

  ( + 0.5  ) * Tai["ei"] * Vijab["nohe"] * Rabcijk["hcdonj"];
  ( - 0.5  ) * Tai["ej"] * Vijab["nohe"] * Rabcijk["hcdoni"];

  ( - 0.5  ) * Tai["dm"] * Vijab["nmgh"] * Rabcijk["hgcnij"];
  ( + 0.5  ) * Tai["cm"] * Vijab["nmgh"] * Rabcijk["hgdnij"];

  ( + 1.0  ) * Tai["en"] * Vijab["onhe"] * Rabcijk["hcdoij"];

  ( + 1.0  ) * Tabcijk["ecdnij"] * Vijab["onhe"] * Rai["ho"];

  ( + 0.5  ) * Tabcijk["efdoij"] * Vijab["poef"] * Rai["cp"];
  ( - 0.5  ) * Tabcijk["efcoij"] * Vijab["poef"] * Rai["dp"];

  ( - 0.5  ) * Tabcijk["ecdnoi"] * Vijab["nohe"] * Rai["hj"];
  ( + 0.5  ) * Tabcijk["ecdnoj"] * Vijab["nohe"] * Rai["hi"];

  //: BEGIN TRIPLES
  ( - 1.0  ) * Fia["oh"] * Tabcijk["hefijk"] * Rai["do"];
  ( - 1.0  ) * Fia["oh"] * Tabcijk["hfdijk"] * Rai["eo"];
  ( - 1.0  ) * Fia["oh"] * Tabcijk["hdeijk"] * Rai["fo"];

  ( - 1.0  ) * Fia["oh"] * Tabcijk["defoij"] * Rai["hk"];
  ( + 1.0  ) * Fia["oh"] * Tabcijk["defoik"] * Rai["hj"];
  ( + 1.0  ) * Fia["oh"] * Tabcijk["defokj"] * Rai["hi"];

  ( - 1.0  ) * Fia["oh"] * Tai["hi"] * Rabcijk["defojk"];
  ( + 1.0  ) * Fia["oh"] * Tai["hj"] * Rabcijk["defoik"];
  ( + 1.0  ) * Fia["oh"] * Tai["hk"] * Rabcijk["defoji"];

  ( - 1.0  ) * Fia["oh"] * Tai["fo"] * Rabcijk["hdeijk"];
  ( + 1.0  ) * Fia["oh"] * Tai["eo"] * Rabcijk["hdfijk"];
  ( + 1.0  ) * Fia["oh"] * Tai["do"] * Rabcijk["hfeijk"];

  ( + 1.0  ) * Fia["oh"] * Tabij["hfij"] * Rabij["deok"];
  ( - 1.0  ) * Fia["oh"] * Tabij["heij"] * Rabij["dfok"];
  ( - 1.0  ) * Fia["oh"] * Tabij["hdij"] * Rabij["feok"];
  ( - 1.0  ) * Fia["oh"] * Tabij["hfik"] * Rabij["deoj"];
  ( + 1.0  ) * Fia["oh"] * Tabij["heik"] * Rabij["dfoj"];
  ( + 1.0  ) * Fia["oh"] * Tabij["hdik"] * Rabij["feoj"];
  ( - 1.0  ) * Fia["oh"] * Tabij["hfkj"] * Rabij["deoi"];
  ( + 1.0  ) * Fia["oh"] * Tabij["hekj"] * Rabij["dfoi"];
  ( + 1.0  ) * Fia["oh"] * Tabij["hdkj"] * Rabij["feoi"];

  ( + 1.0  ) * Fia["oh"] * Tabij["efoi"] * Rabij["hdjk"];
  ( + 1.0  ) * Fia["oh"] * Tabij["fdoi"] * Rabij["hejk"];
  ( + 1.0  ) * Fia["oh"] * Tabij["deoi"] * Rabij["hfjk"];
  ( - 1.0  ) * Fia["oh"] * Tabij["efoj"] * Rabij["hdik"];
  ( - 1.0  ) * Fia["oh"] * Tabij["fdoj"] * Rabij["heik"];
  ( - 1.0  ) * Fia["oh"] * Tabij["deoj"] * Rabij["hfik"];
  ( - 1.0  ) * Fia["oh"] * Tabij["efok"] * Rabij["hdji"];
  ( - 1.0  ) * Fia["oh"] * Tabij["fdok"] * Rabij["heji"];
  ( - 1.0  ) * Fia["oh"] * Tabij["deok"] * Rabij["hfji"];

  ( + 1.0  ) * Viajk["ofij"] * Rabij["deok"];
  ( - 1.0  ) * Viajk["oeij"] * Rabij["dfok"];
  ( - 1.0  ) * Viajk["odij"] * Rabij["feok"];
  ( - 1.0  ) * Viajk["ofik"] * Rabij["deoj"];
  ( + 1.0  ) * Viajk["oeik"] * Rabij["dfoj"];
  ( + 1.0  ) * Viajk["odik"] * Rabij["feoj"];
  ( - 1.0  ) * Viajk["ofkj"] * Rabij["deoi"];
  ( + 1.0  ) * Viajk["oekj"] * Rabij["dfoi"];
  ( + 1.0  ) * Viajk["odkj"] * Rabij["feoi"];

  ( + 1.0  ) * Vabic["efig"] * Rabij["gdjk"];
  ( + 1.0  ) * Vabic["fdig"] * Rabij["gejk"];
  ( + 1.0  ) * Vabic["deig"] * Rabij["gfjk"];
  ( - 1.0  ) * Vabic["efjg"] * Rabij["gdik"];
  ( - 1.0  ) * Vabic["fdjg"] * Rabij["geik"];
  ( - 1.0  ) * Vabic["dejg"] * Rabij["gfik"];
  ( - 1.0  ) * Vabic["efkg"] * Rabij["gdji"];
  ( - 1.0  ) * Vabic["fdkg"] * Rabij["geji"];
  ( - 1.0  ) * Vabic["dekg"] * Rabij["gfji"];

  ( - 0.5  ) * Vijkl["opij"] * Rabcijk["defpok"];
  ( + 0.5  ) * Vijkl["opik"] * Rabcijk["defpoj"];
  ( + 0.5  ) * Vijkl["opkj"] * Rabcijk["defpoi"];

  ( - 1.0  ) * Viajb["ofih"] * Rabcijk["hdeojk"];
  ( + 1.0  ) * Viajb["oeih"] * Rabcijk["hdfojk"];
  ( + 1.0  ) * Viajb["odih"] * Rabcijk["hfeojk"];
  ( + 1.0  ) * Viajb["ofjh"] * Rabcijk["hdeoik"];
  ( - 1.0  ) * Viajb["oejh"] * Rabcijk["hdfoik"];
  ( - 1.0  ) * Viajb["odjh"] * Rabcijk["hfeoik"];
  ( + 1.0  ) * Viajb["ofkh"] * Rabcijk["hdeoji"];
  ( - 1.0  ) * Viajb["oekh"] * Rabcijk["hdfoji"];
  ( - 1.0  ) * Viajb["odkh"] * Rabcijk["hfeoji"];

  ( - 0.5  ) * Vabcd["efgh"] * Rabcijk["hgdijk"];
  ( - 0.5  ) * Vabcd["fdgh"] * Rabcijk["hgeijk"];
  ( - 0.5  ) * Vabcd["degh"] * Rabcijk["hgfijk"];

  ( - 1.0  ) * Tai["fo"] * Vijkl["poij"] * Rabij["depk"];
  ( + 1.0  ) * Tai["eo"] * Vijkl["poij"] * Rabij["dfpk"];
  ( + 1.0  ) * Tai["do"] * Vijkl["poij"] * Rabij["fepk"];
  ( + 1.0  ) * Tai["fo"] * Vijkl["poik"] * Rabij["depj"];
  ( - 1.0  ) * Tai["eo"] * Vijkl["poik"] * Rabij["dfpj"];
  ( - 1.0  ) * Tai["do"] * Vijkl["poik"] * Rabij["fepj"];
  ( + 1.0  ) * Tai["fo"] * Vijkl["pokj"] * Rabij["depi"];
  ( - 1.0  ) * Tai["eo"] * Vijkl["pokj"] * Rabij["dfpi"];
  ( - 1.0  ) * Tai["do"] * Vijkl["pokj"] * Rabij["fepi"];

  ( + 1.0  ) * Tai["gj"] * Viajb["pfig"] * Rabij["depk"];
  ( - 1.0  ) * Tai["gj"] * Viajb["peig"] * Rabij["dfpk"];
  ( - 1.0  ) * Tai["gj"] * Viajb["pdig"] * Rabij["fepk"];
  ( - 1.0  ) * Tai["gk"] * Viajb["pfig"] * Rabij["depj"];
  ( + 1.0  ) * Tai["gk"] * Viajb["peig"] * Rabij["dfpj"];
  ( + 1.0  ) * Tai["gk"] * Viajb["pdig"] * Rabij["fepj"];
  ( - 1.0  ) * Tai["gi"] * Viajb["pfjg"] * Rabij["depk"];
  ( + 1.0  ) * Tai["gi"] * Viajb["pejg"] * Rabij["dfpk"];
  ( + 1.0  ) * Tai["gi"] * Viajb["pdjg"] * Rabij["fepk"];
  ( + 1.0  ) * Tai["gi"] * Viajb["pfkg"] * Rabij["depj"];
  ( - 1.0  ) * Tai["gi"] * Viajb["pekg"] * Rabij["dfpj"];
  ( - 1.0  ) * Tai["gi"] * Viajb["pdkg"] * Rabij["fepj"];
  ( + 1.0  ) * Tai["gk"] * Viajb["pfjg"] * Rabij["depi"];
  ( - 1.0  ) * Tai["gk"] * Viajb["pejg"] * Rabij["dfpi"];
  ( - 1.0  ) * Tai["gk"] * Viajb["pdjg"] * Rabij["fepi"];
  ( - 1.0  ) * Tai["gj"] * Viajb["pfkg"] * Rabij["depi"];
  ( + 1.0  ) * Tai["gj"] * Viajb["pekg"] * Rabij["dfpi"];
  ( + 1.0  ) * Tai["gj"] * Viajb["pdkg"] * Rabij["fepi"];

  ( - 1.0  ) * Tai["eo"] * Viajb["ofih"] * Rabij["hdjk"];
  ( + 1.0  ) * Tai["do"] * Viajb["ofih"] * Rabij["hejk"];
  ( + 1.0  ) * Tai["fo"] * Viajb["oeih"] * Rabij["hdjk"];
  ( - 1.0  ) * Tai["fo"] * Viajb["odih"] * Rabij["hejk"];
  ( - 1.0  ) * Tai["do"] * Viajb["oeih"] * Rabij["hfjk"];
  ( + 1.0  ) * Tai["eo"] * Viajb["odih"] * Rabij["hfjk"];
  ( + 1.0  ) * Tai["eo"] * Viajb["ofjh"] * Rabij["hdik"];
  ( - 1.0  ) * Tai["do"] * Viajb["ofjh"] * Rabij["heik"];
  ( - 1.0  ) * Tai["fo"] * Viajb["oejh"] * Rabij["hdik"];
  ( + 1.0  ) * Tai["fo"] * Viajb["odjh"] * Rabij["heik"];
  ( + 1.0  ) * Tai["do"] * Viajb["oejh"] * Rabij["hfik"];
  ( - 1.0  ) * Tai["eo"] * Viajb["odjh"] * Rabij["hfik"];
  ( + 1.0  ) * Tai["eo"] * Viajb["ofkh"] * Rabij["hdji"];
  ( - 1.0  ) * Tai["do"] * Viajb["ofkh"] * Rabij["heji"];
  ( - 1.0  ) * Tai["fo"] * Viajb["oekh"] * Rabij["hdji"];
  ( + 1.0  ) * Tai["fo"] * Viajb["odkh"] * Rabij["heji"];
  ( + 1.0  ) * Tai["do"] * Viajb["oekh"] * Rabij["hfji"];
  ( - 1.0  ) * Tai["eo"] * Viajb["odkh"] * Rabij["hfji"];

  ( - 1.0  ) * Tai["gi"] * Vabcd["efhg"] * Rabij["hdjk"];
  ( - 1.0  ) * Tai["gi"] * Vabcd["fdhg"] * Rabij["hejk"];
  ( - 1.0  ) * Tai["gi"] * Vabcd["dehg"] * Rabij["hfjk"];
  ( + 1.0  ) * Tai["gj"] * Vabcd["efhg"] * Rabij["hdik"];
  ( + 1.0  ) * Tai["gj"] * Vabcd["fdhg"] * Rabij["heik"];
  ( + 1.0  ) * Tai["gj"] * Vabcd["dehg"] * Rabij["hfik"];
  ( + 1.0  ) * Tai["gk"] * Vabcd["efhg"] * Rabij["hdji"];
  ( + 1.0  ) * Tai["gk"] * Vabcd["fdhg"] * Rabij["heji"];
  ( + 1.0  ) * Tai["gk"] * Vabcd["dehg"] * Rabij["hfji"];

  ( - 0.5  ) * Tai["gj"] * Vijka["pIig"] * Rabcijk["defIpk"];
  ( + 0.5  ) * Tai["gk"] * Vijka["pIig"] * Rabcijk["defIpj"];
  ( + 0.5  ) * Tai["gi"] * Vijka["pIjg"] * Rabcijk["defIpk"];
  ( - 0.5  ) * Tai["gi"] * Vijka["pIkg"] * Rabcijk["defIpj"];
  ( - 0.5  ) * Tai["gk"] * Vijka["pIjg"] * Rabcijk["defIpi"];
  ( + 0.5  ) * Tai["gj"] * Vijka["pIkg"] * Rabcijk["defIpi"];

  ( + 1.0  ) * Tai["fo"] * Vijka["poiA"] * Rabcijk["Adepjk"];
  ( - 1.0  ) * Tai["eo"] * Vijka["poiA"] * Rabcijk["Adfpjk"];
  ( - 1.0  ) * Tai["do"] * Vijka["poiA"] * Rabcijk["Afepjk"];
  ( - 1.0  ) * Tai["fo"] * Vijka["pojA"] * Rabcijk["Adepik"];
  ( + 1.0  ) * Tai["eo"] * Vijka["pojA"] * Rabcijk["Adfpik"];
  ( + 1.0  ) * Tai["do"] * Vijka["pojA"] * Rabcijk["Afepik"];
  ( - 1.0  ) * Tai["fo"] * Vijka["pokA"] * Rabcijk["Adepji"];
  ( + 1.0  ) * Tai["eo"] * Vijka["pokA"] * Rabcijk["Adfpji"];
  ( + 1.0  ) * Tai["do"] * Vijka["pokA"] * Rabcijk["Afepji"];

  ( - 1.0  ) * Tai["gp"] * Vijka["Ipig"] * Rabcijk["defIjk"];
  ( + 1.0  ) * Tai["gp"] * Vijka["Ipjg"] * Rabcijk["defIik"];
  ( + 1.0  ) * Tai["gp"] * Vijka["Ipkg"] * Rabcijk["defIji"];

  ( + 1.0  ) * Tai["gi"] * Viabc["pfAg"] * Rabcijk["Adepjk"];
  ( - 1.0  ) * Tai["gi"] * Viabc["peAg"] * Rabcijk["Adfpjk"];
  ( - 1.0  ) * Tai["gi"] * Viabc["pdAg"] * Rabcijk["Afepjk"];
  ( - 1.0  ) * Tai["gj"] * Viabc["pfAg"] * Rabcijk["Adepik"];
  ( + 1.0  ) * Tai["gj"] * Viabc["peAg"] * Rabcijk["Adfpik"];
  ( + 1.0  ) * Tai["gj"] * Viabc["pdAg"] * Rabcijk["Afepik"];
  ( - 1.0  ) * Tai["gk"] * Viabc["pfAg"] * Rabcijk["Adepji"];
  ( + 1.0  ) * Tai["gk"] * Viabc["peAg"] * Rabcijk["Adfpji"];
  ( + 1.0  ) * Tai["gk"] * Viabc["pdAg"] * Rabcijk["Afepji"];

  ( + 0.5  ) * Tai["eo"] * Viabc["ofhA"] * Rabcijk["Ahdijk"];
  ( - 0.5  ) * Tai["do"] * Viabc["ofhA"] * Rabcijk["Aheijk"];
  ( - 0.5  ) * Tai["fo"] * Viabc["oehA"] * Rabcijk["Ahdijk"];
  ( + 0.5  ) * Tai["fo"] * Viabc["odhA"] * Rabcijk["Aheijk"];
  ( + 0.5  ) * Tai["do"] * Viabc["oehA"] * Rabcijk["Ahfijk"];
  ( - 0.5  ) * Tai["eo"] * Viabc["odhA"] * Rabcijk["Ahfijk"];

  ( - 1.0  ) * Tai["gp"] * Viabc["pfAg"] * Rabcijk["Adeijk"];
  ( + 1.0  ) * Tai["gp"] * Viabc["peAg"] * Rabcijk["Adfijk"];
  ( + 1.0  ) * Tai["gp"] * Viabc["pdAg"] * Rabcijk["Afeijk"];

  ( + 1.0  ) * Tabij["efok"] * Vijkl["poij"] * Rai["dp"];
  ( + 1.0  ) * Tabij["fdok"] * Vijkl["poij"] * Rai["ep"];
  ( + 1.0  ) * Tabij["deok"] * Vijkl["poij"] * Rai["fp"];
  ( - 1.0  ) * Tabij["efoj"] * Vijkl["poik"] * Rai["dp"];
  ( - 1.0  ) * Tabij["fdoj"] * Vijkl["poik"] * Rai["ep"];
  ( - 1.0  ) * Tabij["deoj"] * Vijkl["poik"] * Rai["fp"];
  ( - 1.0  ) * Tabij["efoi"] * Vijkl["pokj"] * Rai["dp"];
  ( - 1.0  ) * Tabij["fdoi"] * Vijkl["pokj"] * Rai["ep"];
  ( - 1.0  ) * Tabij["deoi"] * Vijkl["pokj"] * Rai["fp"];

  ( + 1.0  ) * Tabij["gejk"] * Viajb["pfig"] * Rai["dp"];
  ( - 1.0  ) * Tabij["gdjk"] * Viajb["pfig"] * Rai["ep"];
  ( - 1.0  ) * Tabij["gfjk"] * Viajb["peig"] * Rai["dp"];
  ( + 1.0  ) * Tabij["gfjk"] * Viajb["pdig"] * Rai["ep"];
  ( + 1.0  ) * Tabij["gdjk"] * Viajb["peig"] * Rai["fp"];
  ( - 1.0  ) * Tabij["gejk"] * Viajb["pdig"] * Rai["fp"];
  ( + 1.0  ) * Tabij["geki"] * Viajb["pfjg"] * Rai["dp"];
  ( - 1.0  ) * Tabij["gdki"] * Viajb["pfjg"] * Rai["ep"];
  ( - 1.0  ) * Tabij["gfki"] * Viajb["pejg"] * Rai["dp"];
  ( + 1.0  ) * Tabij["gfki"] * Viajb["pdjg"] * Rai["ep"];
  ( + 1.0  ) * Tabij["gdki"] * Viajb["pejg"] * Rai["fp"];
  ( - 1.0  ) * Tabij["geki"] * Viajb["pdjg"] * Rai["fp"];
  ( + 1.0  ) * Tabij["geij"] * Viajb["pfkg"] * Rai["dp"];
  ( - 1.0  ) * Tabij["gdij"] * Viajb["pfkg"] * Rai["ep"];
  ( - 1.0  ) * Tabij["gfij"] * Viajb["pekg"] * Rai["dp"];
  ( + 1.0  ) * Tabij["gfij"] * Viajb["pdkg"] * Rai["ep"];
  ( + 1.0  ) * Tabij["gdij"] * Viajb["pekg"] * Rai["fp"];
  ( - 1.0  ) * Tabij["geij"] * Viajb["pdkg"] * Rai["fp"];

  ( - 1.0  ) * Tabij["deoj"] * Viajb["ofih"] * Rai["hk"];
  ( + 1.0  ) * Tabij["dfoj"] * Viajb["oeih"] * Rai["hk"];
  ( + 1.0  ) * Tabij["feoj"] * Viajb["odih"] * Rai["hk"];
  ( + 1.0  ) * Tabij["deok"] * Viajb["ofih"] * Rai["hj"];
  ( - 1.0  ) * Tabij["dfok"] * Viajb["oeih"] * Rai["hj"];
  ( - 1.0  ) * Tabij["feok"] * Viajb["odih"] * Rai["hj"];
  ( + 1.0  ) * Tabij["deoi"] * Viajb["ofjh"] * Rai["hk"];
  ( - 1.0  ) * Tabij["dfoi"] * Viajb["oejh"] * Rai["hk"];
  ( - 1.0  ) * Tabij["feoi"] * Viajb["odjh"] * Rai["hk"];
  ( - 1.0  ) * Tabij["deoi"] * Viajb["ofkh"] * Rai["hj"];
  ( + 1.0  ) * Tabij["dfoi"] * Viajb["oekh"] * Rai["hj"];
  ( + 1.0  ) * Tabij["feoi"] * Viajb["odkh"] * Rai["hj"];
  ( - 1.0  ) * Tabij["deok"] * Viajb["ofjh"] * Rai["hi"];
  ( + 1.0  ) * Tabij["dfok"] * Viajb["oejh"] * Rai["hi"];
  ( + 1.0  ) * Tabij["feok"] * Viajb["odjh"] * Rai["hi"];
  ( + 1.0  ) * Tabij["deoj"] * Viajb["ofkh"] * Rai["hi"];
  ( - 1.0  ) * Tabij["dfoj"] * Viajb["oekh"] * Rai["hi"];
  ( - 1.0  ) * Tabij["feoj"] * Viajb["odkh"] * Rai["hi"];

  ( + 1.0  ) * Tabij["gdij"] * Vabcd["efhg"] * Rai["hk"];
  ( - 1.0  ) * Tabij["geij"] * Vabcd["dfhg"] * Rai["hk"];
  ( - 1.0  ) * Tabij["gfij"] * Vabcd["edhg"] * Rai["hk"];
  ( - 1.0  ) * Tabij["gdik"] * Vabcd["efhg"] * Rai["hj"];
  ( + 1.0  ) * Tabij["geik"] * Vabcd["dfhg"] * Rai["hj"];
  ( + 1.0  ) * Tabij["gfik"] * Vabcd["edhg"] * Rai["hj"];
  ( - 1.0  ) * Tabij["gdkj"] * Vabcd["efhg"] * Rai["hi"];
  ( + 1.0  ) * Tabij["gekj"] * Vabcd["dfhg"] * Rai["hi"];
  ( + 1.0  ) * Tabij["gfkj"] * Vabcd["edhg"] * Rai["hi"];

  ( - 0.5  ) * Tabij["gfjk"] * Vijka["pIig"] * Rabij["deIp"];
  ( + 0.5  ) * Tabij["gejk"] * Vijka["pIig"] * Rabij["dfIp"];
  ( + 0.5  ) * Tabij["gdjk"] * Vijka["pIig"] * Rabij["feIp"];
  ( - 0.5  ) * Tabij["gfki"] * Vijka["pIjg"] * Rabij["deIp"];
  ( + 0.5  ) * Tabij["geki"] * Vijka["pIjg"] * Rabij["dfIp"];
  ( + 0.5  ) * Tabij["gdki"] * Vijka["pIjg"] * Rabij["feIp"];
  ( - 0.5  ) * Tabij["gfij"] * Vijka["pIkg"] * Rabij["deIp"];
  ( + 0.5  ) * Tabij["geij"] * Vijka["pIkg"] * Rabij["dfIp"];
  ( + 0.5  ) * Tabij["gdij"] * Vijka["pIkg"] * Rabij["feIp"];

  ( + 1.0  ) * Tabij["efoj"] * Vijka["poiA"] * Rabij["Adpk"];
  ( + 1.0  ) * Tabij["fdoj"] * Vijka["poiA"] * Rabij["Aepk"];
  ( + 1.0  ) * Tabij["deoj"] * Vijka["poiA"] * Rabij["Afpk"];
  ( - 1.0  ) * Tabij["efok"] * Vijka["poiA"] * Rabij["Adpj"];
  ( - 1.0  ) * Tabij["fdok"] * Vijka["poiA"] * Rabij["Aepj"];
  ( - 1.0  ) * Tabij["deok"] * Vijka["poiA"] * Rabij["Afpj"];
  ( - 1.0  ) * Tabij["efoi"] * Vijka["pojA"] * Rabij["Adpk"];
  ( - 1.0  ) * Tabij["fdoi"] * Vijka["pojA"] * Rabij["Aepk"];
  ( - 1.0  ) * Tabij["deoi"] * Vijka["pojA"] * Rabij["Afpk"];
  ( + 1.0  ) * Tabij["efoi"] * Vijka["pokA"] * Rabij["Adpj"];
  ( + 1.0  ) * Tabij["fdoi"] * Vijka["pokA"] * Rabij["Aepj"];
  ( + 1.0  ) * Tabij["deoi"] * Vijka["pokA"] * Rabij["Afpj"];
  ( + 1.0  ) * Tabij["efok"] * Vijka["pojA"] * Rabij["Adpi"];
  ( + 1.0  ) * Tabij["fdok"] * Vijka["pojA"] * Rabij["Aepi"];
  ( + 1.0  ) * Tabij["deok"] * Vijka["pojA"] * Rabij["Afpi"];
  ( - 1.0  ) * Tabij["efoj"] * Vijka["pokA"] * Rabij["Adpi"];
  ( - 1.0  ) * Tabij["fdoj"] * Vijka["pokA"] * Rabij["Aepi"];
  ( - 1.0  ) * Tabij["deoj"] * Vijka["pokA"] * Rabij["Afpi"];

  ( + 1.0  ) * Tabij["gfpj"] * Vijka["Ipig"] * Rabij["deIk"];
  ( - 1.0  ) * Tabij["gepj"] * Vijka["Ipig"] * Rabij["dfIk"];
  ( - 1.0  ) * Tabij["gdpj"] * Vijka["Ipig"] * Rabij["feIk"];
  ( - 1.0  ) * Tabij["gfpk"] * Vijka["Ipig"] * Rabij["deIj"];
  ( + 1.0  ) * Tabij["gepk"] * Vijka["Ipig"] * Rabij["dfIj"];
  ( + 1.0  ) * Tabij["gdpk"] * Vijka["Ipig"] * Rabij["feIj"];
  ( - 1.0  ) * Tabij["gfpi"] * Vijka["Ipjg"] * Rabij["deIk"];
  ( + 1.0  ) * Tabij["gepi"] * Vijka["Ipjg"] * Rabij["dfIk"];
  ( + 1.0  ) * Tabij["gdpi"] * Vijka["Ipjg"] * Rabij["feIk"];
  ( + 1.0  ) * Tabij["gfpi"] * Vijka["Ipkg"] * Rabij["deIj"];
  ( - 1.0  ) * Tabij["gepi"] * Vijka["Ipkg"] * Rabij["dfIj"];
  ( - 1.0  ) * Tabij["gdpi"] * Vijka["Ipkg"] * Rabij["feIj"];
  ( + 1.0  ) * Tabij["gfpk"] * Vijka["Ipjg"] * Rabij["deIi"];
  ( - 1.0  ) * Tabij["gepk"] * Vijka["Ipjg"] * Rabij["dfIi"];
  ( - 1.0  ) * Tabij["gdpk"] * Vijka["Ipjg"] * Rabij["feIi"];
  ( - 1.0  ) * Tabij["gfpj"] * Vijka["Ipkg"] * Rabij["deIi"];
  ( + 1.0  ) * Tabij["gepj"] * Vijka["Ipkg"] * Rabij["dfIi"];
  ( + 1.0  ) * Tabij["gdpj"] * Vijka["Ipkg"] * Rabij["feIi"];

  ( + 0.5  ) * Tabij["efop"] * Vijka["opiA"] * Rabij["Adjk"];
  ( + 0.5  ) * Tabij["fdop"] * Vijka["opiA"] * Rabij["Aejk"];
  ( + 0.5  ) * Tabij["deop"] * Vijka["opiA"] * Rabij["Afjk"];
  ( - 0.5  ) * Tabij["efop"] * Vijka["opjA"] * Rabij["Adik"];
  ( - 0.5  ) * Tabij["fdop"] * Vijka["opjA"] * Rabij["Aeik"];
  ( - 0.5  ) * Tabij["deop"] * Vijka["opjA"] * Rabij["Afik"];
  ( - 0.5  ) * Tabij["efop"] * Vijka["opkA"] * Rabij["Adji"];
  ( - 0.5  ) * Tabij["fdop"] * Vijka["opkA"] * Rabij["Aeji"];
  ( - 0.5  ) * Tabij["deop"] * Vijka["opkA"] * Rabij["Afji"];

  ( - 1.0  ) * Tabij["geij"] * Viabc["pfAg"] * Rabij["Adpk"];
  ( + 1.0  ) * Tabij["gdij"] * Viabc["pfAg"] * Rabij["Aepk"];
  ( + 1.0  ) * Tabij["gfij"] * Viabc["peAg"] * Rabij["Adpk"];
  ( - 1.0  ) * Tabij["gfij"] * Viabc["pdAg"] * Rabij["Aepk"];
  ( - 1.0  ) * Tabij["gdij"] * Viabc["peAg"] * Rabij["Afpk"];
  ( + 1.0  ) * Tabij["geij"] * Viabc["pdAg"] * Rabij["Afpk"];
  ( + 1.0  ) * Tabij["geik"] * Viabc["pfAg"] * Rabij["Adpj"];
  ( - 1.0  ) * Tabij["gdik"] * Viabc["pfAg"] * Rabij["Aepj"];
  ( - 1.0  ) * Tabij["gfik"] * Viabc["peAg"] * Rabij["Adpj"];
  ( + 1.0  ) * Tabij["gfik"] * Viabc["pdAg"] * Rabij["Aepj"];
  ( + 1.0  ) * Tabij["gdik"] * Viabc["peAg"] * Rabij["Afpj"];
  ( - 1.0  ) * Tabij["geik"] * Viabc["pdAg"] * Rabij["Afpj"];
  ( + 1.0  ) * Tabij["gekj"] * Viabc["pfAg"] * Rabij["Adpi"];
  ( - 1.0  ) * Tabij["gdkj"] * Viabc["pfAg"] * Rabij["Aepi"];
  ( - 1.0  ) * Tabij["gfkj"] * Viabc["peAg"] * Rabij["Adpi"];
  ( + 1.0  ) * Tabij["gfkj"] * Viabc["pdAg"] * Rabij["Aepi"];
  ( + 1.0  ) * Tabij["gdkj"] * Viabc["peAg"] * Rabij["Afpi"];
  ( - 1.0  ) * Tabij["gekj"] * Viabc["pdAg"] * Rabij["Afpi"];

  ( + 0.5  ) * Tabij["ghij"] * Viabc["Ifgh"] * Rabij["deIk"];
  ( - 0.5  ) * Tabij["ghij"] * Viabc["Iegh"] * Rabij["dfIk"];
  ( - 0.5  ) * Tabij["ghij"] * Viabc["Idgh"] * Rabij["feIk"];
  ( - 0.5  ) * Tabij["ghik"] * Viabc["Ifgh"] * Rabij["deIj"];
  ( + 0.5  ) * Tabij["ghik"] * Viabc["Iegh"] * Rabij["dfIj"];
  ( + 0.5  ) * Tabij["ghik"] * Viabc["Idgh"] * Rabij["feIj"];
  ( - 0.5  ) * Tabij["ghkj"] * Viabc["Ifgh"] * Rabij["deIi"];
  ( + 0.5  ) * Tabij["ghkj"] * Viabc["Iegh"] * Rabij["dfIi"];
  ( + 0.5  ) * Tabij["ghkj"] * Viabc["Idgh"] * Rabij["feIi"];

  ( - 0.5  ) * Tabij["deoi"] * Viabc["ofhA"] * Rabij["Ahjk"];
  ( + 0.5  ) * Tabij["dfoi"] * Viabc["oehA"] * Rabij["Ahjk"];
  ( + 0.5  ) * Tabij["feoi"] * Viabc["odhA"] * Rabij["Ahjk"];
  ( + 0.5  ) * Tabij["deoj"] * Viabc["ofhA"] * Rabij["Ahik"];
  ( - 0.5  ) * Tabij["dfoj"] * Viabc["oehA"] * Rabij["Ahik"];
  ( - 0.5  ) * Tabij["feoj"] * Viabc["odhA"] * Rabij["Ahik"];
  ( + 0.5  ) * Tabij["deok"] * Viabc["ofhA"] * Rabij["Ahji"];
  ( - 0.5  ) * Tabij["dfok"] * Viabc["oehA"] * Rabij["Ahji"];
  ( - 0.5  ) * Tabij["feok"] * Viabc["odhA"] * Rabij["Ahji"];

  ( - 1.0  ) * Tabij["gepi"] * Viabc["pfAg"] * Rabij["Adjk"];
  ( + 1.0  ) * Tabij["gdpi"] * Viabc["pfAg"] * Rabij["Aejk"];
  ( + 1.0  ) * Tabij["gfpi"] * Viabc["peAg"] * Rabij["Adjk"];
  ( - 1.0  ) * Tabij["gfpi"] * Viabc["pdAg"] * Rabij["Aejk"];
  ( - 1.0  ) * Tabij["gdpi"] * Viabc["peAg"] * Rabij["Afjk"];
  ( + 1.0  ) * Tabij["gepi"] * Viabc["pdAg"] * Rabij["Afjk"];
  ( + 1.0  ) * Tabij["gepj"] * Viabc["pfAg"] * Rabij["Adik"];
  ( - 1.0  ) * Tabij["gdpj"] * Viabc["pfAg"] * Rabij["Aeik"];
  ( - 1.0  ) * Tabij["gfpj"] * Viabc["peAg"] * Rabij["Adik"];
  ( + 1.0  ) * Tabij["gfpj"] * Viabc["pdAg"] * Rabij["Aeik"];
  ( + 1.0  ) * Tabij["gdpj"] * Viabc["peAg"] * Rabij["Afik"];
  ( - 1.0  ) * Tabij["gepj"] * Viabc["pdAg"] * Rabij["Afik"];
  ( + 1.0  ) * Tabij["gepk"] * Viabc["pfAg"] * Rabij["Adji"];
  ( - 1.0  ) * Tabij["gdpk"] * Viabc["pfAg"] * Rabij["Aeji"];
  ( - 1.0  ) * Tabij["gfpk"] * Viabc["peAg"] * Rabij["Adji"];
  ( + 1.0  ) * Tabij["gfpk"] * Viabc["pdAg"] * Rabij["Aeji"];
  ( + 1.0  ) * Tabij["gdpk"] * Viabc["peAg"] * Rabij["Afji"];
  ( - 1.0  ) * Tabij["gepk"] * Viabc["pdAg"] * Rabij["Afji"];

  ( - 0.5  ) * Tabij["gfij"] * Vijab["pIBg"] * Rabcijk["BdeIpk"];
  ( + 0.5  ) * Tabij["geij"] * Vijab["pIBg"] * Rabcijk["BdfIpk"];
  ( + 0.5  ) * Tabij["gdij"] * Vijab["pIBg"] * Rabcijk["BfeIpk"];
  ( + 0.5  ) * Tabij["gfik"] * Vijab["pIBg"] * Rabcijk["BdeIpj"];
  ( - 0.5  ) * Tabij["geik"] * Vijab["pIBg"] * Rabcijk["BdfIpj"];
  ( - 0.5  ) * Tabij["gdik"] * Vijab["pIBg"] * Rabcijk["BfeIpj"];
  ( + 0.5  ) * Tabij["gfkj"] * Vijab["pIBg"] * Rabcijk["BdeIpi"];
  ( - 0.5  ) * Tabij["gekj"] * Vijab["pIBg"] * Rabcijk["BdfIpi"];
  ( - 0.5  ) * Tabij["gdkj"] * Vijab["pIBg"] * Rabcijk["BfeIpi"];

  ( - 0.25  ) * Tabij["ghij"] * Vijab["IJgh"] * Rabcijk["defJIk"];
  ( + 0.25  ) * Tabij["ghik"] * Vijab["IJgh"] * Rabcijk["defJIj"];
  ( + 0.25  ) * Tabij["ghkj"] * Vijab["IJgh"] * Rabcijk["defJIi"];

  ( - 0.5  ) * Tabij["efoi"] * Vijab["poAB"] * Rabcijk["BAdpjk"];
  ( - 0.5  ) * Tabij["fdoi"] * Vijab["poAB"] * Rabcijk["BAepjk"];
  ( - 0.5  ) * Tabij["deoi"] * Vijab["poAB"] * Rabcijk["BAfpjk"];
  ( + 0.5  ) * Tabij["efoj"] * Vijab["poAB"] * Rabcijk["BAdpik"];
  ( + 0.5  ) * Tabij["fdoj"] * Vijab["poAB"] * Rabcijk["BAepik"];
  ( + 0.5  ) * Tabij["deoj"] * Vijab["poAB"] * Rabcijk["BAfpik"];
  ( + 0.5  ) * Tabij["efok"] * Vijab["poAB"] * Rabcijk["BAdpji"];
  ( + 0.5  ) * Tabij["fdok"] * Vijab["poAB"] * Rabcijk["BAepji"];
  ( + 0.5  ) * Tabij["deok"] * Vijab["poAB"] * Rabcijk["BAfpji"];

  ( + 1.0  ) * Tabij["gfpi"] * Vijab["IpBg"] * Rabcijk["BdeIjk"];
  ( - 1.0  ) * Tabij["gepi"] * Vijab["IpBg"] * Rabcijk["BdfIjk"];
  ( - 1.0  ) * Tabij["gdpi"] * Vijab["IpBg"] * Rabcijk["BfeIjk"];
  ( - 1.0  ) * Tabij["gfpj"] * Vijab["IpBg"] * Rabcijk["BdeIik"];
  ( + 1.0  ) * Tabij["gepj"] * Vijab["IpBg"] * Rabcijk["BdfIik"];
  ( + 1.0  ) * Tabij["gdpj"] * Vijab["IpBg"] * Rabcijk["BfeIik"];
  ( - 1.0  ) * Tabij["gfpk"] * Vijab["IpBg"] * Rabcijk["BdeIji"];
  ( + 1.0  ) * Tabij["gepk"] * Vijab["IpBg"] * Rabcijk["BdfIji"];
  ( + 1.0  ) * Tabij["gdpk"] * Vijab["IpBg"] * Rabcijk["BfeIji"];

  ( + 0.5  ) * Tabij["ghIi"] * Vijab["JIgh"] * Rabcijk["defJjk"];
  ( - 0.5  ) * Tabij["ghIj"] * Vijab["JIgh"] * Rabcijk["defJik"];
  ( - 0.5  ) * Tabij["ghIk"] * Vijab["JIgh"] * Rabcijk["defJji"];

  ( - 0.25  ) * Tabij["efop"] * Vijab["opAB"] * Rabcijk["BAdijk"];
  ( - 0.25  ) * Tabij["fdop"] * Vijab["opAB"] * Rabcijk["BAeijk"];
  ( - 0.25  ) * Tabij["deop"] * Vijab["opAB"] * Rabcijk["BAfijk"];

  ( + 0.5  ) * Tabij["gfpI"] * Vijab["pIBg"] * Rabcijk["Bdeijk"];
  ( - 0.5  ) * Tabij["gepI"] * Vijab["pIBg"] * Rabcijk["Bdfijk"];
  ( - 0.5  ) * Tabij["gdpI"] * Vijab["pIBg"] * Rabcijk["Bfeijk"];

  ( + 1.0  ) * Tabcijk["defojk"] * Vijka["poiA"] * Rai["Ap"];
  ( + 1.0  ) * Tabcijk["defoki"] * Vijka["pojA"] * Rai["Ap"];
  ( + 1.0  ) * Tabcijk["defoij"] * Vijka["pokA"] * Rai["Ap"];

  ( - 1.0  ) * Tabcijk["gefpjk"] * Vijka["Ipig"] * Rai["dI"];
  ( - 1.0  ) * Tabcijk["gfdpjk"] * Vijka["Ipig"] * Rai["eI"];
  ( - 1.0  ) * Tabcijk["gdepjk"] * Vijka["Ipig"] * Rai["fI"];
  ( - 1.0  ) * Tabcijk["gefpki"] * Vijka["Ipjg"] * Rai["dI"];
  ( - 1.0  ) * Tabcijk["gfdpki"] * Vijka["Ipjg"] * Rai["eI"];
  ( - 1.0  ) * Tabcijk["gdepki"] * Vijka["Ipjg"] * Rai["fI"];
  ( - 1.0  ) * Tabcijk["gefpij"] * Vijka["Ipkg"] * Rai["dI"];
  ( - 1.0  ) * Tabcijk["gfdpij"] * Vijka["Ipkg"] * Rai["eI"];
  ( - 1.0  ) * Tabcijk["gdepij"] * Vijka["Ipkg"] * Rai["fI"];

  ( - 0.5  ) * Tabcijk["defopj"] * Vijka["opiA"] * Rai["Ak"];
  ( + 0.5  ) * Tabcijk["defopk"] * Vijka["opiA"] * Rai["Aj"];
  ( + 0.5  ) * Tabcijk["defopi"] * Vijka["opjA"] * Rai["Ak"];
  ( - 0.5  ) * Tabcijk["defopi"] * Vijka["opkA"] * Rai["Aj"];
  ( - 0.5  ) * Tabcijk["defopk"] * Vijka["opjA"] * Rai["Ai"];
  ( + 0.5  ) * Tabcijk["defopj"] * Vijka["opkA"] * Rai["Ai"];

  ( + 1.0  ) * Tabcijk["gdeijk"] * Viabc["pfAg"] * Rai["Ap"];
  ( - 1.0  ) * Tabcijk["gdfijk"] * Viabc["peAg"] * Rai["Ap"];
  ( - 1.0  ) * Tabcijk["gfeijk"] * Viabc["pdAg"] * Rai["Ap"];

  ( + 0.5  ) * Tabcijk["gheijk"] * Viabc["Ifgh"] * Rai["dI"];
  ( - 0.5  ) * Tabcijk["ghdijk"] * Viabc["Ifgh"] * Rai["eI"];
  ( - 0.5  ) * Tabcijk["ghfijk"] * Viabc["Iegh"] * Rai["dI"];
  ( + 0.5  ) * Tabcijk["ghfijk"] * Viabc["Idgh"] * Rai["eI"];
  ( + 0.5  ) * Tabcijk["ghdijk"] * Viabc["Iegh"] * Rai["fI"];
  ( - 0.5  ) * Tabcijk["gheijk"] * Viabc["Idgh"] * Rai["fI"];

  ( - 1.0  ) * Tabcijk["gdepij"] * Viabc["pfAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabcijk["gdfpij"] * Viabc["peAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabcijk["gfepij"] * Viabc["pdAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabcijk["gdepik"] * Viabc["pfAg"] * Rai["Aj"];
  ( - 1.0  ) * Tabcijk["gdfpik"] * Viabc["peAg"] * Rai["Aj"];
  ( - 1.0  ) * Tabcijk["gfepik"] * Viabc["pdAg"] * Rai["Aj"];
  ( + 1.0  ) * Tabcijk["gdepkj"] * Viabc["pfAg"] * Rai["Ai"];
  ( - 1.0  ) * Tabcijk["gdfpkj"] * Viabc["peAg"] * Rai["Ai"];
  ( - 1.0  ) * Tabcijk["gfepkj"] * Viabc["pdAg"] * Rai["Ai"];

  ( + 0.5  ) * Tabcijk["gefijk"] * Vijab["pIBg"] * Rabij["BdIp"];
  ( + 0.5  ) * Tabcijk["gfdijk"] * Vijab["pIBg"] * Rabij["BeIp"];
  ( + 0.5  ) * Tabcijk["gdeijk"] * Vijab["pIBg"] * Rabij["BfIp"];

  ( - 0.25  ) * Tabcijk["ghfijk"] * Vijab["IJgh"] * Rabij["deJI"];
  ( + 0.25  ) * Tabcijk["gheijk"] * Vijab["IJgh"] * Rabij["dfJI"];
  ( + 0.25  ) * Tabcijk["ghdijk"] * Vijab["IJgh"] * Rabij["feJI"];

  ( + 0.5  ) * Tabcijk["defoij"] * Vijab["poAB"] * Rabij["BApk"];
  ( - 0.5  ) * Tabcijk["defoik"] * Vijab["poAB"] * Rabij["BApj"];
  ( - 0.5  ) * Tabcijk["defokj"] * Vijab["poAB"] * Rabij["BApi"];

  ( + 1.0  ) * Tabcijk["gefpij"] * Vijab["IpBg"] * Rabij["BdIk"];
  ( + 1.0  ) * Tabcijk["gfdpij"] * Vijab["IpBg"] * Rabij["BeIk"];
  ( + 1.0  ) * Tabcijk["gdepij"] * Vijab["IpBg"] * Rabij["BfIk"];
  ( - 1.0  ) * Tabcijk["gefpik"] * Vijab["IpBg"] * Rabij["BdIj"];
  ( - 1.0  ) * Tabcijk["gfdpik"] * Vijab["IpBg"] * Rabij["BeIj"];
  ( - 1.0  ) * Tabcijk["gdepik"] * Vijab["IpBg"] * Rabij["BfIj"];
  ( - 1.0  ) * Tabcijk["gefpkj"] * Vijab["IpBg"] * Rabij["BdIi"];
  ( - 1.0  ) * Tabcijk["gfdpkj"] * Vijab["IpBg"] * Rabij["BeIi"];
  ( - 1.0  ) * Tabcijk["gdepkj"] * Vijab["IpBg"] * Rabij["BfIi"];

  ( - 0.5  ) * Tabcijk["ghfIij"] * Vijab["JIgh"] * Rabij["deJk"];
  ( + 0.5  ) * Tabcijk["gheIij"] * Vijab["JIgh"] * Rabij["dfJk"];
  ( + 0.5  ) * Tabcijk["ghdIij"] * Vijab["JIgh"] * Rabij["feJk"];
  ( + 0.5  ) * Tabcijk["ghfIik"] * Vijab["JIgh"] * Rabij["deJj"];
  ( - 0.5  ) * Tabcijk["gheIik"] * Vijab["JIgh"] * Rabij["dfJj"];
  ( - 0.5  ) * Tabcijk["ghdIik"] * Vijab["JIgh"] * Rabij["feJj"];
  ( + 0.5  ) * Tabcijk["ghfIkj"] * Vijab["JIgh"] * Rabij["deJi"];
  ( - 0.5  ) * Tabcijk["gheIkj"] * Vijab["JIgh"] * Rabij["dfJi"];
  ( - 0.5  ) * Tabcijk["ghdIkj"] * Vijab["JIgh"] * Rabij["feJi"];

  ( - 0.25  ) * Tabcijk["defopi"] * Vijab["opAB"] * Rabij["BAjk"];
  ( + 0.25  ) * Tabcijk["defopj"] * Vijab["opAB"] * Rabij["BAik"];
  ( + 0.25  ) * Tabcijk["defopk"] * Vijab["opAB"] * Rabij["BAji"];

  ( - 0.5  ) * Tabcijk["gefpIi"] * Vijab["pIBg"] * Rabij["Bdjk"];
  ( - 0.5  ) * Tabcijk["gfdpIi"] * Vijab["pIBg"] * Rabij["Bejk"];
  ( - 0.5  ) * Tabcijk["gdepIi"] * Vijab["pIBg"] * Rabij["Bfjk"];
  ( + 0.5  ) * Tabcijk["gefpIj"] * Vijab["pIBg"] * Rabij["Bdik"];
  ( + 0.5  ) * Tabcijk["gfdpIj"] * Vijab["pIBg"] * Rabij["Beik"];
  ( + 0.5  ) * Tabcijk["gdepIj"] * Vijab["pIBg"] * Rabij["Bfik"];
  ( + 0.5  ) * Tabcijk["gefpIk"] * Vijab["pIBg"] * Rabij["Bdji"];
  ( + 0.5  ) * Tabcijk["gfdpIk"] * Vijab["pIBg"] * Rabij["Beji"];
  ( + 0.5  ) * Tabcijk["gdepIk"] * Vijab["pIBg"] * Rabij["Bfji"];

  ( - 1.0  ) * Tai["fo"] * Tai["hj"] * Vijka["Ioih"] * Rabij["deIk"];
  ( + 1.0  ) * Tai["eo"] * Tai["hj"] * Vijka["Ioih"] * Rabij["dfIk"];
  ( + 1.0  ) * Tai["do"] * Tai["hj"] * Vijka["Ioih"] * Rabij["feIk"];
  ( + 1.0  ) * Tai["fo"] * Tai["hk"] * Vijka["Ioih"] * Rabij["deIj"];
  ( - 1.0  ) * Tai["eo"] * Tai["hk"] * Vijka["Ioih"] * Rabij["dfIj"];
  ( - 1.0  ) * Tai["do"] * Tai["hk"] * Vijka["Ioih"] * Rabij["feIj"];
  ( + 1.0  ) * Tai["fo"] * Tai["hi"] * Vijka["Iojh"] * Rabij["deIk"];
  ( - 1.0  ) * Tai["eo"] * Tai["hi"] * Vijka["Iojh"] * Rabij["dfIk"];
  ( - 1.0  ) * Tai["do"] * Tai["hi"] * Vijka["Iojh"] * Rabij["feIk"];
  ( - 1.0  ) * Tai["fo"] * Tai["hi"] * Vijka["Iokh"] * Rabij["deIj"];
  ( + 1.0  ) * Tai["eo"] * Tai["hi"] * Vijka["Iokh"] * Rabij["dfIj"];
  ( + 1.0  ) * Tai["do"] * Tai["hi"] * Vijka["Iokh"] * Rabij["feIj"];
  ( - 1.0  ) * Tai["fo"] * Tai["hk"] * Vijka["Iojh"] * Rabij["deIi"];
  ( + 1.0  ) * Tai["eo"] * Tai["hk"] * Vijka["Iojh"] * Rabij["dfIi"];
  ( + 1.0  ) * Tai["do"] * Tai["hk"] * Vijka["Iojh"] * Rabij["feIi"];
  ( + 1.0  ) * Tai["fo"] * Tai["hj"] * Vijka["Iokh"] * Rabij["deIi"];
  ( - 1.0  ) * Tai["eo"] * Tai["hj"] * Vijka["Iokh"] * Rabij["dfIi"];
  ( - 1.0  ) * Tai["do"] * Tai["hj"] * Vijka["Iokh"] * Rabij["feIi"];

  ( + 0.5  ) * Tai["fo"] * Tai["ep"] * Vijka["poiA"] * Rabij["Adjk"];
  ( - 0.5  ) * Tai["eo"] * Tai["fp"] * Vijka["poiA"] * Rabij["Adjk"];
  ( - 0.5  ) * Tai["fo"] * Tai["dp"] * Vijka["poiA"] * Rabij["Aejk"];
  ( + 0.5  ) * Tai["do"] * Tai["fp"] * Vijka["poiA"] * Rabij["Aejk"];
  ( + 0.5  ) * Tai["eo"] * Tai["dp"] * Vijka["poiA"] * Rabij["Afjk"];
  ( - 0.5  ) * Tai["do"] * Tai["ep"] * Vijka["poiA"] * Rabij["Afjk"];
  ( - 0.5  ) * Tai["fo"] * Tai["ep"] * Vijka["pojA"] * Rabij["Adik"];
  ( + 0.5  ) * Tai["eo"] * Tai["fp"] * Vijka["pojA"] * Rabij["Adik"];
  ( + 0.5  ) * Tai["fo"] * Tai["dp"] * Vijka["pojA"] * Rabij["Aeik"];
  ( - 0.5  ) * Tai["do"] * Tai["fp"] * Vijka["pojA"] * Rabij["Aeik"];
  ( - 0.5  ) * Tai["eo"] * Tai["dp"] * Vijka["pojA"] * Rabij["Afik"];
  ( + 0.5  ) * Tai["do"] * Tai["ep"] * Vijka["pojA"] * Rabij["Afik"];
  ( - 0.5  ) * Tai["fo"] * Tai["ep"] * Vijka["pokA"] * Rabij["Adji"];
  ( + 0.5  ) * Tai["eo"] * Tai["fp"] * Vijka["pokA"] * Rabij["Adji"];
  ( + 0.5  ) * Tai["fo"] * Tai["dp"] * Vijka["pokA"] * Rabij["Aeji"];
  ( - 0.5  ) * Tai["do"] * Tai["fp"] * Vijka["pokA"] * Rabij["Aeji"];
  ( - 0.5  ) * Tai["eo"] * Tai["dp"] * Vijka["pokA"] * Rabij["Afji"];
  ( + 0.5  ) * Tai["do"] * Tai["ep"] * Vijka["pokA"] * Rabij["Afji"];

  ( - 0.5  ) * Tai["gi"] * Tai["hj"] * Viabc["Ifhg"] * Rabij["deIk"];
  ( + 0.5  ) * Tai["gj"] * Tai["hi"] * Viabc["Ifhg"] * Rabij["deIk"];
  ( + 0.5  ) * Tai["gi"] * Tai["hj"] * Viabc["Iehg"] * Rabij["dfIk"];
  ( - 0.5  ) * Tai["gj"] * Tai["hi"] * Viabc["Iehg"] * Rabij["dfIk"];
  ( + 0.5  ) * Tai["gi"] * Tai["hj"] * Viabc["Idhg"] * Rabij["feIk"];
  ( - 0.5  ) * Tai["gj"] * Tai["hi"] * Viabc["Idhg"] * Rabij["feIk"];
  ( + 0.5  ) * Tai["gi"] * Tai["hk"] * Viabc["Ifhg"] * Rabij["deIj"];
  ( - 0.5  ) * Tai["gk"] * Tai["hi"] * Viabc["Ifhg"] * Rabij["deIj"];
  ( - 0.5  ) * Tai["gi"] * Tai["hk"] * Viabc["Iehg"] * Rabij["dfIj"];
  ( + 0.5  ) * Tai["gk"] * Tai["hi"] * Viabc["Iehg"] * Rabij["dfIj"];
  ( - 0.5  ) * Tai["gi"] * Tai["hk"] * Viabc["Idhg"] * Rabij["feIj"];
  ( + 0.5  ) * Tai["gk"] * Tai["hi"] * Viabc["Idhg"] * Rabij["feIj"];
  ( - 0.5  ) * Tai["gj"] * Tai["hk"] * Viabc["Ifhg"] * Rabij["deIi"];
  ( + 0.5  ) * Tai["gk"] * Tai["hj"] * Viabc["Ifhg"] * Rabij["deIi"];
  ( + 0.5  ) * Tai["gj"] * Tai["hk"] * Viabc["Iehg"] * Rabij["dfIi"];
  ( - 0.5  ) * Tai["gk"] * Tai["hj"] * Viabc["Iehg"] * Rabij["dfIi"];
  ( + 0.5  ) * Tai["gj"] * Tai["hk"] * Viabc["Idhg"] * Rabij["feIi"];
  ( - 0.5  ) * Tai["gk"] * Tai["hj"] * Viabc["Idhg"] * Rabij["feIi"];

  ( + 1.0  ) * Tai["eo"] * Tai["hi"] * Viabc["ofAh"] * Rabij["Adjk"];
  ( - 1.0  ) * Tai["do"] * Tai["hi"] * Viabc["ofAh"] * Rabij["Aejk"];
  ( - 1.0  ) * Tai["fo"] * Tai["hi"] * Viabc["oeAh"] * Rabij["Adjk"];
  ( + 1.0  ) * Tai["fo"] * Tai["hi"] * Viabc["odAh"] * Rabij["Aejk"];
  ( + 1.0  ) * Tai["do"] * Tai["hi"] * Viabc["oeAh"] * Rabij["Afjk"];
  ( - 1.0  ) * Tai["eo"] * Tai["hi"] * Viabc["odAh"] * Rabij["Afjk"];
  ( - 1.0  ) * Tai["eo"] * Tai["hj"] * Viabc["ofAh"] * Rabij["Adik"];
  ( + 1.0  ) * Tai["do"] * Tai["hj"] * Viabc["ofAh"] * Rabij["Aeik"];
  ( + 1.0  ) * Tai["fo"] * Tai["hj"] * Viabc["oeAh"] * Rabij["Adik"];
  ( - 1.0  ) * Tai["fo"] * Tai["hj"] * Viabc["odAh"] * Rabij["Aeik"];
  ( - 1.0  ) * Tai["do"] * Tai["hj"] * Viabc["oeAh"] * Rabij["Afik"];
  ( + 1.0  ) * Tai["eo"] * Tai["hj"] * Viabc["odAh"] * Rabij["Afik"];
  ( - 1.0  ) * Tai["eo"] * Tai["hk"] * Viabc["ofAh"] * Rabij["Adji"];
  ( + 1.0  ) * Tai["do"] * Tai["hk"] * Viabc["ofAh"] * Rabij["Aeji"];
  ( + 1.0  ) * Tai["fo"] * Tai["hk"] * Viabc["oeAh"] * Rabij["Adji"];
  ( - 1.0  ) * Tai["fo"] * Tai["hk"] * Viabc["odAh"] * Rabij["Aeji"];
  ( - 1.0  ) * Tai["do"] * Tai["hk"] * Viabc["oeAh"] * Rabij["Afji"];
  ( + 1.0  ) * Tai["eo"] * Tai["hk"] * Viabc["odAh"] * Rabij["Afji"];

  ( + 0.25  ) * Tai["gi"] * Tai["hj"] * Vijab["IJhg"] * Rabcijk["defJIk"];
  ( - 0.25  ) * Tai["gj"] * Tai["hi"] * Vijab["IJhg"] * Rabcijk["defJIk"];
  ( - 0.25  ) * Tai["gi"] * Tai["hk"] * Vijab["IJhg"] * Rabcijk["defJIj"];
  ( + 0.25  ) * Tai["gk"] * Tai["hi"] * Vijab["IJhg"] * Rabcijk["defJIj"];
  ( + 0.25  ) * Tai["gj"] * Tai["hk"] * Vijab["IJhg"] * Rabcijk["defJIi"];
  ( - 0.25  ) * Tai["gk"] * Tai["hj"] * Vijab["IJhg"] * Rabcijk["defJIi"];

  ( - 1.0  ) * Tai["fo"] * Tai["hi"] * Vijab["IoBh"] * Rabcijk["BdeIjk"];
  ( + 1.0  ) * Tai["eo"] * Tai["hi"] * Vijab["IoBh"] * Rabcijk["BdfIjk"];
  ( + 1.0  ) * Tai["do"] * Tai["hi"] * Vijab["IoBh"] * Rabcijk["BfeIjk"];
  ( + 1.0  ) * Tai["fo"] * Tai["hj"] * Vijab["IoBh"] * Rabcijk["BdeIik"];
  ( - 1.0  ) * Tai["eo"] * Tai["hj"] * Vijab["IoBh"] * Rabcijk["BdfIik"];
  ( - 1.0  ) * Tai["do"] * Tai["hj"] * Vijab["IoBh"] * Rabcijk["BfeIik"];
  ( + 1.0  ) * Tai["fo"] * Tai["hk"] * Vijab["IoBh"] * Rabcijk["BdeIji"];
  ( - 1.0  ) * Tai["eo"] * Tai["hk"] * Vijab["IoBh"] * Rabcijk["BdfIji"];
  ( - 1.0  ) * Tai["do"] * Tai["hk"] * Vijab["IoBh"] * Rabcijk["BfeIji"];

  ( + 1.0  ) * Tai["gi"] * Tai["hI"] * Vijab["JIhg"] * Rabcijk["defJjk"];
  ( - 1.0  ) * Tai["gj"] * Tai["hI"] * Vijab["JIhg"] * Rabcijk["defJik"];
  ( - 1.0  ) * Tai["gk"] * Tai["hI"] * Vijab["JIhg"] * Rabcijk["defJji"];

  ( - 0.25  ) * Tai["fo"] * Tai["ep"] * Vijab["poAB"] * Rabcijk["BAdijk"];
  ( + 0.25  ) * Tai["eo"] * Tai["fp"] * Vijab["poAB"] * Rabcijk["BAdijk"];
  ( + 0.25  ) * Tai["fo"] * Tai["dp"] * Vijab["poAB"] * Rabcijk["BAeijk"];
  ( - 0.25  ) * Tai["do"] * Tai["fp"] * Vijab["poAB"] * Rabcijk["BAeijk"];
  ( - 0.25  ) * Tai["eo"] * Tai["dp"] * Vijab["poAB"] * Rabcijk["BAfijk"];
  ( + 0.25  ) * Tai["do"] * Tai["ep"] * Vijab["poAB"] * Rabcijk["BAfijk"];

  ( + 1.0  ) * Tai["fo"] * Tai["hI"] * Vijab["IoBh"] * Rabcijk["Bdeijk"];
  ( - 1.0  ) * Tai["eo"] * Tai["hI"] * Vijab["IoBh"] * Rabcijk["Bdfijk"];
  ( - 1.0  ) * Tai["do"] * Tai["hI"] * Vijab["IoBh"] * Rabcijk["Bfeijk"];

  ( + 1.0  ) * Tabij["efok"] * Tai["hj"] * Vijka["Ioih"] * Rai["dI"];
  ( + 1.0  ) * Tabij["fdok"] * Tai["hj"] * Vijka["Ioih"] * Rai["eI"];
  ( + 1.0  ) * Tabij["deok"] * Tai["hj"] * Vijka["Ioih"] * Rai["fI"];
  ( - 1.0  ) * Tabij["efoj"] * Tai["hk"] * Vijka["Ioih"] * Rai["dI"];
  ( - 1.0  ) * Tabij["fdoj"] * Tai["hk"] * Vijka["Ioih"] * Rai["eI"];
  ( - 1.0  ) * Tabij["deoj"] * Tai["hk"] * Vijka["Ioih"] * Rai["fI"];
  ( - 1.0  ) * Tabij["efok"] * Tai["hi"] * Vijka["Iojh"] * Rai["dI"];
  ( - 1.0  ) * Tabij["fdok"] * Tai["hi"] * Vijka["Iojh"] * Rai["eI"];
  ( - 1.0  ) * Tabij["deok"] * Tai["hi"] * Vijka["Iojh"] * Rai["fI"];
  ( + 1.0  ) * Tabij["efoj"] * Tai["hi"] * Vijka["Iokh"] * Rai["dI"];
  ( + 1.0  ) * Tabij["fdoj"] * Tai["hi"] * Vijka["Iokh"] * Rai["eI"];
  ( + 1.0  ) * Tabij["deoj"] * Tai["hi"] * Vijka["Iokh"] * Rai["fI"];
  ( + 1.0  ) * Tabij["efoi"] * Tai["hk"] * Vijka["Iojh"] * Rai["dI"];
  ( + 1.0  ) * Tabij["fdoi"] * Tai["hk"] * Vijka["Iojh"] * Rai["eI"];
  ( + 1.0  ) * Tabij["deoi"] * Tai["hk"] * Vijka["Iojh"] * Rai["fI"];
  ( - 1.0  ) * Tabij["efoi"] * Tai["hj"] * Vijka["Iokh"] * Rai["dI"];
  ( - 1.0  ) * Tabij["fdoi"] * Tai["hj"] * Vijka["Iokh"] * Rai["eI"];
  ( - 1.0  ) * Tabij["deoi"] * Tai["hj"] * Vijka["Iokh"] * Rai["fI"];

  ( - 1.0  ) * Tabij["gejk"] * Tai["fp"] * Vijka["Ipig"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gdjk"] * Tai["fp"] * Vijka["Ipig"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gfjk"] * Tai["ep"] * Vijka["Ipig"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gfjk"] * Tai["dp"] * Vijka["Ipig"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gdjk"] * Tai["ep"] * Vijka["Ipig"] * Rai["fI"];
  ( + 1.0  ) * Tabij["gejk"] * Tai["dp"] * Vijka["Ipig"] * Rai["fI"];
  ( - 1.0  ) * Tabij["geki"] * Tai["fp"] * Vijka["Ipjg"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gdki"] * Tai["fp"] * Vijka["Ipjg"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gfki"] * Tai["ep"] * Vijka["Ipjg"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gfki"] * Tai["dp"] * Vijka["Ipjg"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gdki"] * Tai["ep"] * Vijka["Ipjg"] * Rai["fI"];
  ( + 1.0  ) * Tabij["geki"] * Tai["dp"] * Vijka["Ipjg"] * Rai["fI"];
  ( - 1.0  ) * Tabij["geij"] * Tai["fp"] * Vijka["Ipkg"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gdij"] * Tai["fp"] * Vijka["Ipkg"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gfij"] * Tai["ep"] * Vijka["Ipkg"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gfij"] * Tai["dp"] * Vijka["Ipkg"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gdij"] * Tai["ep"] * Vijka["Ipkg"] * Rai["fI"];
  ( + 1.0  ) * Tabij["geij"] * Tai["dp"] * Vijka["Ipkg"] * Rai["fI"];

  ( - 1.0  ) * Tabij["deoj"] * Tai["fp"] * Vijka["poiA"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["dfoj"] * Tai["ep"] * Vijka["poiA"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["feoj"] * Tai["dp"] * Vijka["poiA"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["deok"] * Tai["fp"] * Vijka["poiA"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["dfok"] * Tai["ep"] * Vijka["poiA"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["feok"] * Tai["dp"] * Vijka["poiA"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["deoi"] * Tai["fp"] * Vijka["pojA"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["dfoi"] * Tai["ep"] * Vijka["pojA"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["feoi"] * Tai["dp"] * Vijka["pojA"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["deoi"] * Tai["fp"] * Vijka["pokA"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["dfoi"] * Tai["ep"] * Vijka["pokA"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["feoi"] * Tai["dp"] * Vijka["pokA"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["deok"] * Tai["fp"] * Vijka["pojA"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["dfok"] * Tai["ep"] * Vijka["pojA"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["feok"] * Tai["dp"] * Vijka["pojA"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["deoj"] * Tai["fp"] * Vijka["pokA"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["dfoj"] * Tai["ep"] * Vijka["pokA"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["feoj"] * Tai["dp"] * Vijka["pokA"] * Rai["Ai"];

  ( + 1.0  ) * Tabij["gejk"] * Tai["hi"] * Viabc["Ifhg"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gdjk"] * Tai["hi"] * Viabc["Ifhg"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gfjk"] * Tai["hi"] * Viabc["Iehg"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gfjk"] * Tai["hi"] * Viabc["Idhg"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gdjk"] * Tai["hi"] * Viabc["Iehg"] * Rai["fI"];
  ( - 1.0  ) * Tabij["gejk"] * Tai["hi"] * Viabc["Idhg"] * Rai["fI"];
  ( + 1.0  ) * Tabij["geki"] * Tai["hj"] * Viabc["Ifhg"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gdki"] * Tai["hj"] * Viabc["Ifhg"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gfki"] * Tai["hj"] * Viabc["Iehg"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gfki"] * Tai["hj"] * Viabc["Idhg"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gdki"] * Tai["hj"] * Viabc["Iehg"] * Rai["fI"];
  ( - 1.0  ) * Tabij["geki"] * Tai["hj"] * Viabc["Idhg"] * Rai["fI"];
  ( + 1.0  ) * Tabij["geij"] * Tai["hk"] * Viabc["Ifhg"] * Rai["dI"];
  ( - 1.0  ) * Tabij["gdij"] * Tai["hk"] * Viabc["Ifhg"] * Rai["eI"];
  ( - 1.0  ) * Tabij["gfij"] * Tai["hk"] * Viabc["Iehg"] * Rai["dI"];
  ( + 1.0  ) * Tabij["gfij"] * Tai["hk"] * Viabc["Idhg"] * Rai["eI"];
  ( + 1.0  ) * Tabij["gdij"] * Tai["hk"] * Viabc["Iehg"] * Rai["fI"];
  ( - 1.0  ) * Tabij["geij"] * Tai["hk"] * Viabc["Idhg"] * Rai["fI"];

  ( + 1.0  ) * Tabij["deoj"] * Tai["hi"] * Viabc["ofAh"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["dfoj"] * Tai["hi"] * Viabc["oeAh"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["feoj"] * Tai["hi"] * Viabc["odAh"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["deok"] * Tai["hi"] * Viabc["ofAh"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["dfok"] * Tai["hi"] * Viabc["oeAh"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["feok"] * Tai["hi"] * Viabc["odAh"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["deoi"] * Tai["hj"] * Viabc["ofAh"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["dfoi"] * Tai["hj"] * Viabc["oeAh"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["feoi"] * Tai["hj"] * Viabc["odAh"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["deoi"] * Tai["hk"] * Viabc["ofAh"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["dfoi"] * Tai["hk"] * Viabc["oeAh"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["feoi"] * Tai["hk"] * Viabc["odAh"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["deok"] * Tai["hj"] * Viabc["ofAh"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["dfok"] * Tai["hj"] * Viabc["oeAh"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["feok"] * Tai["hj"] * Viabc["odAh"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["deoj"] * Tai["hk"] * Viabc["ofAh"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["dfoj"] * Tai["hk"] * Viabc["oeAh"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["feoj"] * Tai["hk"] * Viabc["odAh"] * Rai["Ai"];

  ( - 1.0  ) * Tabij["gdij"] * Tai["ep"] * Viabc["pfAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["geij"] * Tai["dp"] * Viabc["pfAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["gdij"] * Tai["fp"] * Viabc["peAg"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["geij"] * Tai["fp"] * Viabc["pdAg"] * Rai["Ak"];
  ( - 1.0  ) * Tabij["gfij"] * Tai["dp"] * Viabc["peAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["gfij"] * Tai["ep"] * Viabc["pdAg"] * Rai["Ak"];
  ( + 1.0  ) * Tabij["gdik"] * Tai["ep"] * Viabc["pfAg"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["geik"] * Tai["dp"] * Viabc["pfAg"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["gdik"] * Tai["fp"] * Viabc["peAg"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["geik"] * Tai["fp"] * Viabc["pdAg"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["gfik"] * Tai["dp"] * Viabc["peAg"] * Rai["Aj"];
  ( - 1.0  ) * Tabij["gfik"] * Tai["ep"] * Viabc["pdAg"] * Rai["Aj"];
  ( + 1.0  ) * Tabij["gdkj"] * Tai["ep"] * Viabc["pfAg"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["gekj"] * Tai["dp"] * Viabc["pfAg"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["gdkj"] * Tai["fp"] * Viabc["peAg"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["gekj"] * Tai["fp"] * Viabc["pdAg"] * Rai["Ai"];
  ( + 1.0  ) * Tabij["gfkj"] * Tai["dp"] * Viabc["peAg"] * Rai["Ai"];
  ( - 1.0  ) * Tabij["gfkj"] * Tai["ep"] * Viabc["pdAg"] * Rai["Ai"];

  ( - 0.5  ) * Tabij["gfjk"] * Tai["hi"] * Vijab["IJhg"] * Rabij["deJI"];
  ( + 0.5  ) * Tabij["gejk"] * Tai["hi"] * Vijab["IJhg"] * Rabij["dfJI"];
  ( + 0.5  ) * Tabij["gdjk"] * Tai["hi"] * Vijab["IJhg"] * Rabij["feJI"];
  ( - 0.5  ) * Tabij["gfki"] * Tai["hj"] * Vijab["IJhg"] * Rabij["deJI"];
  ( + 0.5  ) * Tabij["geki"] * Tai["hj"] * Vijab["IJhg"] * Rabij["dfJI"];
  ( + 0.5  ) * Tabij["gdki"] * Tai["hj"] * Vijab["IJhg"] * Rabij["feJI"];
  ( - 0.5  ) * Tabij["gfij"] * Tai["hk"] * Vijab["IJhg"] * Rabij["deJI"];
  ( + 0.5  ) * Tabij["geij"] * Tai["hk"] * Vijab["IJhg"] * Rabij["dfJI"];
  ( + 0.5  ) * Tabij["gdij"] * Tai["hk"] * Vijab["IJhg"] * Rabij["feJI"];

  ( - 1.0  ) * Tabij["efoj"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BdIk"];
  ( - 1.0  ) * Tabij["fdoj"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BeIk"];
  ( - 1.0  ) * Tabij["deoj"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BfIk"];
  ( + 1.0  ) * Tabij["efok"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BdIj"];
  ( + 1.0  ) * Tabij["fdok"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BeIj"];
  ( + 1.0  ) * Tabij["deok"] * Tai["hi"] * Vijab["IoBh"] * Rabij["BfIj"];
  ( + 1.0  ) * Tabij["efoi"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BdIk"];
  ( + 1.0  ) * Tabij["fdoi"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BeIk"];
  ( + 1.0  ) * Tabij["deoi"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BfIk"];
  ( - 1.0  ) * Tabij["efoi"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BdIj"];
  ( - 1.0  ) * Tabij["fdoi"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BeIj"];
  ( - 1.0  ) * Tabij["deoi"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BfIj"];
  ( - 1.0  ) * Tabij["efok"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BdIi"];
  ( - 1.0  ) * Tabij["fdok"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BeIi"];
  ( - 1.0  ) * Tabij["deok"] * Tai["hj"] * Vijab["IoBh"] * Rabij["BfIi"];
  ( + 1.0  ) * Tabij["efoj"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BdIi"];
  ( + 1.0  ) * Tabij["fdoj"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BeIi"];
  ( + 1.0  ) * Tabij["deoj"] * Tai["hk"] * Vijab["IoBh"] * Rabij["BfIi"];

  ( + 1.0  ) * Tabij["gfpj"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["deJk"];
  ( - 1.0  ) * Tabij["gepj"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["dfJk"];
  ( - 1.0  ) * Tabij["gdpj"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["feJk"];
  ( - 1.0  ) * Tabij["gfpk"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["deJj"];
  ( + 1.0  ) * Tabij["gepk"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["dfJj"];
  ( + 1.0  ) * Tabij["gdpk"] * Tai["Ai"] * Vijab["JpAg"] * Rabij["feJj"];
  ( - 1.0  ) * Tabij["gfpi"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["deJk"];
  ( + 1.0  ) * Tabij["gepi"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["dfJk"];
  ( + 1.0  ) * Tabij["gdpi"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["feJk"];
  ( + 1.0  ) * Tabij["gfpi"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["deJj"];
  ( - 1.0  ) * Tabij["gepi"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["dfJj"];
  ( - 1.0  ) * Tabij["gdpi"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["feJj"];
  ( + 1.0  ) * Tabij["gfpk"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["deJi"];
  ( - 1.0  ) * Tabij["gepk"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["dfJi"];
  ( - 1.0  ) * Tabij["gdpk"] * Tai["Aj"] * Vijab["JpAg"] * Rabij["feJi"];
  ( - 1.0  ) * Tabij["gfpj"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["deJi"];
  ( + 1.0  ) * Tabij["gepj"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["dfJi"];
  ( + 1.0  ) * Tabij["gdpj"] * Tai["Ak"] * Vijab["JpAg"] * Rabij["feJi"];

  ( - 0.5  ) * Tabij["efop"] * Tai["Ai"] * Vijab["opBA"] * Rabij["Bdjk"];
  ( - 0.5  ) * Tabij["fdop"] * Tai["Ai"] * Vijab["opBA"] * Rabij["Bejk"];
  ( - 0.5  ) * Tabij["deop"] * Tai["Ai"] * Vijab["opBA"] * Rabij["Bfjk"];
  ( + 0.5  ) * Tabij["efop"] * Tai["Aj"] * Vijab["opBA"] * Rabij["Bdik"];
  ( + 0.5  ) * Tabij["fdop"] * Tai["Aj"] * Vijab["opBA"] * Rabij["Beik"];
  ( + 0.5  ) * Tabij["deop"] * Tai["Aj"] * Vijab["opBA"] * Rabij["Bfik"];
  ( + 0.5  ) * Tabij["efop"] * Tai["Ak"] * Vijab["opBA"] * Rabij["Bdji"];
  ( + 0.5  ) * Tabij["fdop"] * Tai["Ak"] * Vijab["opBA"] * Rabij["Beji"];
  ( + 0.5  ) * Tabij["deop"] * Tai["Ak"] * Vijab["opBA"] * Rabij["Bfji"];

  ( + 1.0  ) * Tabij["geij"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BdIk"];
  ( - 1.0  ) * Tabij["gdij"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BeIk"];
  ( - 1.0  ) * Tabij["gfij"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BdIk"];
  ( + 1.0  ) * Tabij["gfij"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BeIk"];
  ( + 1.0  ) * Tabij["gdij"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BfIk"];
  ( - 1.0  ) * Tabij["geij"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BfIk"];
  ( - 1.0  ) * Tabij["geik"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BdIj"];
  ( + 1.0  ) * Tabij["gdik"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BeIj"];
  ( + 1.0  ) * Tabij["gfik"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BdIj"];
  ( - 1.0  ) * Tabij["gfik"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BeIj"];
  ( - 1.0  ) * Tabij["gdik"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BfIj"];
  ( + 1.0  ) * Tabij["geik"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BfIj"];
  ( - 1.0  ) * Tabij["gekj"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BdIi"];
  ( + 1.0  ) * Tabij["gdkj"] * Tai["fp"] * Vijab["IpBg"] * Rabij["BeIi"];
  ( + 1.0  ) * Tabij["gfkj"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BdIi"];
  ( - 1.0  ) * Tabij["gfkj"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BeIi"];
  ( - 1.0  ) * Tabij["gdkj"] * Tai["ep"] * Vijab["IpBg"] * Rabij["BfIi"];
  ( + 1.0  ) * Tabij["gekj"] * Tai["dp"] * Vijab["IpBg"] * Rabij["BfIi"];

  ( - 0.5  ) * Tabij["ghij"] * Tai["fI"] * Vijab["JIgh"] * Rabij["deJk"];
  ( + 0.5  ) * Tabij["ghij"] * Tai["eI"] * Vijab["JIgh"] * Rabij["dfJk"];
  ( + 0.5  ) * Tabij["ghij"] * Tai["dI"] * Vijab["JIgh"] * Rabij["feJk"];
  ( + 0.5  ) * Tabij["ghik"] * Tai["fI"] * Vijab["JIgh"] * Rabij["deJj"];
  ( - 0.5  ) * Tabij["ghik"] * Tai["eI"] * Vijab["JIgh"] * Rabij["dfJj"];
  ( - 0.5  ) * Tabij["ghik"] * Tai["dI"] * Vijab["JIgh"] * Rabij["feJj"];
  ( + 0.5  ) * Tabij["ghkj"] * Tai["fI"] * Vijab["JIgh"] * Rabij["deJi"];
  ( - 0.5  ) * Tabij["ghkj"] * Tai["eI"] * Vijab["JIgh"] * Rabij["dfJi"];
  ( - 0.5  ) * Tabij["ghkj"] * Tai["dI"] * Vijab["JIgh"] * Rabij["feJi"];

  ( - 1.0  ) * Tabij["gfij"] * Tai["hI"] * Vijab["JIhg"] * Rabij["deJk"];
  ( + 1.0  ) * Tabij["geij"] * Tai["hI"] * Vijab["JIhg"] * Rabij["dfJk"];
  ( + 1.0  ) * Tabij["gdij"] * Tai["hI"] * Vijab["JIhg"] * Rabij["feJk"];
  ( + 1.0  ) * Tabij["gfik"] * Tai["hI"] * Vijab["JIhg"] * Rabij["deJj"];
  ( - 1.0  ) * Tabij["geik"] * Tai["hI"] * Vijab["JIhg"] * Rabij["dfJj"];
  ( - 1.0  ) * Tabij["gdik"] * Tai["hI"] * Vijab["JIhg"] * Rabij["feJj"];
  ( + 1.0  ) * Tabij["gfkj"] * Tai["hI"] * Vijab["JIhg"] * Rabij["deJi"];
  ( - 1.0  ) * Tabij["gekj"] * Tai["hI"] * Vijab["JIhg"] * Rabij["dfJi"];
  ( - 1.0  ) * Tabij["gdkj"] * Tai["hI"] * Vijab["JIhg"] * Rabij["feJi"];

  ( - 0.5  ) * Tabij["deoi"] * Tai["fp"] * Vijab["poAB"] * Rabij["BAjk"];
  ( + 0.5  ) * Tabij["dfoi"] * Tai["ep"] * Vijab["poAB"] * Rabij["BAjk"];
  ( + 0.5  ) * Tabij["feoi"] * Tai["dp"] * Vijab["poAB"] * Rabij["BAjk"];
  ( + 0.5  ) * Tabij["deoj"] * Tai["fp"] * Vijab["poAB"] * Rabij["BAik"];
  ( - 0.5  ) * Tabij["dfoj"] * Tai["ep"] * Vijab["poAB"] * Rabij["BAik"];
  ( - 0.5  ) * Tabij["feoj"] * Tai["dp"] * Vijab["poAB"] * Rabij["BAik"];
  ( + 0.5  ) * Tabij["deok"] * Tai["fp"] * Vijab["poAB"] * Rabij["BAji"];
  ( - 0.5  ) * Tabij["dfok"] * Tai["ep"] * Vijab["poAB"] * Rabij["BAji"];
  ( - 0.5  ) * Tabij["feok"] * Tai["dp"] * Vijab["poAB"] * Rabij["BAji"];

  ( - 1.0  ) * Tabij["gepi"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Bdjk"];
  ( + 1.0  ) * Tabij["gdpi"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Bejk"];
  ( + 1.0  ) * Tabij["gfpi"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bdjk"];
  ( - 1.0  ) * Tabij["gfpi"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Bejk"];
  ( - 1.0  ) * Tabij["gdpi"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bfjk"];
  ( + 1.0  ) * Tabij["gepi"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Bfjk"];
  ( + 1.0  ) * Tabij["gepj"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Bdik"];
  ( - 1.0  ) * Tabij["gdpj"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Beik"];
  ( - 1.0  ) * Tabij["gfpj"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bdik"];
  ( + 1.0  ) * Tabij["gfpj"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Beik"];
  ( + 1.0  ) * Tabij["gdpj"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bfik"];
  ( - 1.0  ) * Tabij["gepj"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Bfik"];
  ( + 1.0  ) * Tabij["gepk"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Bdji"];
  ( - 1.0  ) * Tabij["gdpk"] * Tai["fI"] * Vijab["IpBg"] * Rabij["Beji"];
  ( - 1.0  ) * Tabij["gfpk"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bdji"];
  ( + 1.0  ) * Tabij["gfpk"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Beji"];
  ( + 1.0  ) * Tabij["gdpk"] * Tai["eI"] * Vijab["IpBg"] * Rabij["Bfji"];
  ( - 1.0  ) * Tabij["gepk"] * Tai["dI"] * Vijab["IpBg"] * Rabij["Bfji"];

  ( - 1.0  ) * Tabij["efoi"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bdjk"];
  ( - 1.0  ) * Tabij["fdoi"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bejk"];
  ( - 1.0  ) * Tabij["deoi"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bfjk"];
  ( + 1.0  ) * Tabij["efoj"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bdik"];
  ( + 1.0  ) * Tabij["fdoj"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Beik"];
  ( + 1.0  ) * Tabij["deoj"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bfik"];
  ( + 1.0  ) * Tabij["efok"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bdji"];
  ( + 1.0  ) * Tabij["fdok"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Beji"];
  ( + 1.0  ) * Tabij["deok"] * Tai["hI"] * Vijab["IoBh"] * Rabij["Bfji"];

  ( - 1.0  ) * Tabcijk["defojk"] * Tai["hi"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabcijk["defoki"] * Tai["hj"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabcijk["defoij"] * Tai["hk"] * Vijab["IoBh"] * Rai["BI"];

  ( - 1.0  ) * Tabcijk["gefpjk"] * Tai["Ai"] * Vijab["JpAg"] * Rai["dJ"];
  ( - 1.0  ) * Tabcijk["gfdpjk"] * Tai["Ai"] * Vijab["JpAg"] * Rai["eJ"];
  ( - 1.0  ) * Tabcijk["gdepjk"] * Tai["Ai"] * Vijab["JpAg"] * Rai["fJ"];
  ( - 1.0  ) * Tabcijk["gefpki"] * Tai["Aj"] * Vijab["JpAg"] * Rai["dJ"];
  ( - 1.0  ) * Tabcijk["gfdpki"] * Tai["Aj"] * Vijab["JpAg"] * Rai["eJ"];
  ( - 1.0  ) * Tabcijk["gdepki"] * Tai["Aj"] * Vijab["JpAg"] * Rai["fJ"];
  ( - 1.0  ) * Tabcijk["gefpij"] * Tai["Ak"] * Vijab["JpAg"] * Rai["dJ"];
  ( - 1.0  ) * Tabcijk["gfdpij"] * Tai["Ak"] * Vijab["JpAg"] * Rai["eJ"];
  ( - 1.0  ) * Tabcijk["gdepij"] * Tai["Ak"] * Vijab["JpAg"] * Rai["fJ"];

  ( + 0.5  ) * Tabcijk["defopj"] * Tai["Ai"] * Vijab["opBA"] * Rai["Bk"];
  ( - 0.5  ) * Tabcijk["defopk"] * Tai["Ai"] * Vijab["opBA"] * Rai["Bj"];
  ( - 0.5  ) * Tabcijk["defopi"] * Tai["Aj"] * Vijab["opBA"] * Rai["Bk"];
  ( + 0.5  ) * Tabcijk["defopi"] * Tai["Ak"] * Vijab["opBA"] * Rai["Bj"];
  ( + 0.5  ) * Tabcijk["defopk"] * Tai["Aj"] * Vijab["opBA"] * Rai["Bi"];
  ( - 0.5  ) * Tabcijk["defopj"] * Tai["Ak"] * Vijab["opBA"] * Rai["Bi"];

  ( - 1.0  ) * Tabcijk["gdeijk"] * Tai["fp"] * Vijab["IpBg"] * Rai["BI"];
  ( + 1.0  ) * Tabcijk["gdfijk"] * Tai["ep"] * Vijab["IpBg"] * Rai["BI"];
  ( + 1.0  ) * Tabcijk["gfeijk"] * Tai["dp"] * Vijab["IpBg"] * Rai["BI"];

  ( - 0.5  ) * Tabcijk["gheijk"] * Tai["fI"] * Vijab["JIgh"] * Rai["dJ"];
  ( + 0.5  ) * Tabcijk["ghdijk"] * Tai["fI"] * Vijab["JIgh"] * Rai["eJ"];
  ( + 0.5  ) * Tabcijk["ghfijk"] * Tai["eI"] * Vijab["JIgh"] * Rai["dJ"];
  ( - 0.5  ) * Tabcijk["ghfijk"] * Tai["dI"] * Vijab["JIgh"] * Rai["eJ"];
  ( - 0.5  ) * Tabcijk["ghdijk"] * Tai["eI"] * Vijab["JIgh"] * Rai["fJ"];
  ( + 0.5  ) * Tabcijk["gheijk"] * Tai["dI"] * Vijab["JIgh"] * Rai["fJ"];

  ( + 1.0  ) * Tabcijk["gefijk"] * Tai["hI"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 1.0  ) * Tabcijk["gfdijk"] * Tai["hI"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 1.0  ) * Tabcijk["gdeijk"] * Tai["hI"] * Vijab["JIhg"] * Rai["fJ"];

  ( - 1.0  ) * Tabcijk["gdepij"] * Tai["fI"] * Vijab["IpBg"] * Rai["Bk"];
  ( + 1.0  ) * Tabcijk["gdfpij"] * Tai["eI"] * Vijab["IpBg"] * Rai["Bk"];
  ( + 1.0  ) * Tabcijk["gfepij"] * Tai["dI"] * Vijab["IpBg"] * Rai["Bk"];
  ( + 1.0  ) * Tabcijk["gdepik"] * Tai["fI"] * Vijab["IpBg"] * Rai["Bj"];
  ( - 1.0  ) * Tabcijk["gdfpik"] * Tai["eI"] * Vijab["IpBg"] * Rai["Bj"];
  ( - 1.0  ) * Tabcijk["gfepik"] * Tai["dI"] * Vijab["IpBg"] * Rai["Bj"];
  ( + 1.0  ) * Tabcijk["gdepkj"] * Tai["fI"] * Vijab["IpBg"] * Rai["Bi"];
  ( - 1.0  ) * Tabcijk["gdfpkj"] * Tai["eI"] * Vijab["IpBg"] * Rai["Bi"];
  ( - 1.0  ) * Tabcijk["gfepkj"] * Tai["dI"] * Vijab["IpBg"] * Rai["Bi"];

  ( + 1.0  ) * Tabcijk["defoij"] * Tai["hI"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tabcijk["defoik"] * Tai["hI"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tabcijk["defokj"] * Tai["hI"] * Vijab["IoBh"] * Rai["Bi"];

  ( + 1.0  ) * Tabij["deok"] * Tabij["hfij"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabij["dfok"] * Tabij["heij"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabij["feok"] * Tabij["hdij"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabij["deoj"] * Tabij["hfik"] * Vijab["IoBh"] * Rai["BI"];
  ( + 1.0  ) * Tabij["dfoj"] * Tabij["heik"] * Vijab["IoBh"] * Rai["BI"];
  ( + 1.0  ) * Tabij["feoj"] * Tabij["hdik"] * Vijab["IoBh"] * Rai["BI"];
  ( + 1.0  ) * Tabij["feoi"] * Tabij["hdkj"] * Vijab["IoBh"] * Rai["BI"];
  ( + 1.0  ) * Tabij["dfoi"] * Tabij["hekj"] * Vijab["IoBh"] * Rai["BI"];
  ( - 1.0  ) * Tabij["deoi"] * Tabij["hfkj"] * Vijab["IoBh"] * Rai["BI"];

  ( - 1.0  ) * Tabij["gfij"] * Tabij["heIk"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 1.0  ) * Tabij["gfij"] * Tabij["hdIk"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 1.0  ) * Tabij["geij"] * Tabij["hfIk"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 1.0  ) * Tabij["gdij"] * Tabij["hfIk"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 1.0  ) * Tabij["geij"] * Tabij["hdIk"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 1.0  ) * Tabij["gdij"] * Tabij["heIk"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 1.0  ) * Tabij["gfik"] * Tabij["heIj"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 1.0  ) * Tabij["gfik"] * Tabij["hdIj"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 1.0  ) * Tabij["geik"] * Tabij["hfIj"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 1.0  ) * Tabij["gdik"] * Tabij["hfIj"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 1.0  ) * Tabij["geik"] * Tabij["hdIj"] * Vijab["JIhg"] * Rai["fJ"];
  ( - 1.0  ) * Tabij["gdik"] * Tabij["heIj"] * Vijab["JIhg"] * Rai["fJ"];
  ( - 1.0  ) * Tabij["gekj"] * Tabij["hfIi"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 1.0  ) * Tabij["gdkj"] * Tabij["hfIi"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 1.0  ) * Tabij["gfkj"] * Tabij["heIi"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 1.0  ) * Tabij["gfkj"] * Tabij["hdIi"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 1.0  ) * Tabij["gdkj"] * Tabij["heIi"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 1.0  ) * Tabij["gekj"] * Tabij["hdIi"] * Vijab["JIhg"] * Rai["fJ"];

  ( + 0.5  ) * Tabij["efok"] * Tabij["hAij"] * Vijab["JohA"] * Rai["dJ"];
  ( + 0.5  ) * Tabij["fdok"] * Tabij["hAij"] * Vijab["JohA"] * Rai["eJ"];
  ( + 0.5  ) * Tabij["deok"] * Tabij["hAij"] * Vijab["JohA"] * Rai["fJ"];
  ( - 0.5  ) * Tabij["efoj"] * Tabij["hAik"] * Vijab["JohA"] * Rai["dJ"];
  ( - 0.5  ) * Tabij["fdoj"] * Tabij["hAik"] * Vijab["JohA"] * Rai["eJ"];
  ( - 0.5  ) * Tabij["deoj"] * Tabij["hAik"] * Vijab["JohA"] * Rai["fJ"];
  ( - 0.5  ) * Tabij["efoi"] * Tabij["hAkj"] * Vijab["JohA"] * Rai["dJ"];
  ( - 0.5  ) * Tabij["fdoi"] * Tabij["hAkj"] * Vijab["JohA"] * Rai["eJ"];
  ( - 0.5  ) * Tabij["deoi"] * Tabij["hAkj"] * Vijab["JohA"] * Rai["fJ"];

  ( + 0.5  ) * Tabij["gfij"] * Tabij["depI"] * Vijab["pIBg"] * Rai["Bk"];
  ( - 0.5  ) * Tabij["geij"] * Tabij["dfpI"] * Vijab["pIBg"] * Rai["Bk"];
  ( - 0.5  ) * Tabij["gdij"] * Tabij["fepI"] * Vijab["pIBg"] * Rai["Bk"];
  ( - 0.5  ) * Tabij["gfik"] * Tabij["depI"] * Vijab["pIBg"] * Rai["Bj"];
  ( + 0.5  ) * Tabij["geik"] * Tabij["dfpI"] * Vijab["pIBg"] * Rai["Bj"];
  ( + 0.5  ) * Tabij["gdik"] * Tabij["fepI"] * Vijab["pIBg"] * Rai["Bj"];
  ( - 0.5  ) * Tabij["gfkj"] * Tabij["depI"] * Vijab["pIBg"] * Rai["Bi"];
  ( + 0.5  ) * Tabij["gekj"] * Tabij["dfpI"] * Vijab["pIBg"] * Rai["Bi"];
  ( + 0.5  ) * Tabij["gdkj"] * Tabij["fepI"] * Vijab["pIBg"] * Rai["Bi"];

  ( + 1.0  ) * Tabij["efoi"] * Tabij["hdIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( + 1.0  ) * Tabij["fdoi"] * Tabij["heIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tabij["deoj"] * Tabij["hfIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( + 1.0  ) * Tabij["deoi"] * Tabij["hfIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tabij["fdoj"] * Tabij["heIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tabij["efoj"] * Tabij["hdIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tabij["efoi"] * Tabij["hdIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tabij["fdoi"] * Tabij["heIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tabij["deok"] * Tabij["hfIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tabij["deoi"] * Tabij["hfIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tabij["fdok"] * Tabij["heIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tabij["efok"] * Tabij["hdIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tabij["efoj"] * Tabij["hdIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( + 1.0  ) * Tabij["fdoj"] * Tabij["heIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( - 1.0  ) * Tabij["deok"] * Tabij["hfIj"] * Vijab["IoBh"] * Rai["Bi"];
  ( + 1.0  ) * Tabij["deoj"] * Tabij["hfIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( - 1.0  ) * Tabij["fdok"] * Tabij["heIj"] * Vijab["IoBh"] * Rai["Bi"];
  ( - 1.0  ) * Tabij["efok"] * Tabij["hdIj"] * Vijab["IoBh"] * Rai["Bi"];

  ( + 0.5  ) * Tai["fo"] * Tai["hi"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["deJk"];
  ( - 0.5  ) * Tai["fo"] * Tai["hj"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["deJk"];
  ( - 0.5  ) * Tai["eo"] * Tai["hi"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["dfJk"];
  ( + 0.5  ) * Tai["eo"] * Tai["hj"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["dfJk"];
  ( - 0.5  ) * Tai["do"] * Tai["hi"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["feJk"];
  ( + 0.5  ) * Tai["do"] * Tai["hj"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["feJk"];
  ( - 0.5  ) * Tai["fo"] * Tai["hi"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["deJj"];
  ( + 0.5  ) * Tai["fo"] * Tai["hk"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["deJj"];
  ( + 0.5  ) * Tai["eo"] * Tai["hi"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["dfJj"];
  ( - 0.5  ) * Tai["eo"] * Tai["hk"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["dfJj"];
  ( + 0.5  ) * Tai["do"] * Tai["hi"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["feJj"];
  ( - 0.5  ) * Tai["do"] * Tai["hk"] * Tai["Ai"] * Vijab["JoAh"] * Rabij["feJj"];
  ( + 0.5  ) * Tai["fo"] * Tai["hj"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["deJi"];
  ( - 0.5  ) * Tai["fo"] * Tai["hk"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["deJi"];
  ( - 0.5  ) * Tai["eo"] * Tai["hj"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["dfJi"];
  ( + 0.5  ) * Tai["eo"] * Tai["hk"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["dfJi"];
  ( - 0.5  ) * Tai["do"] * Tai["hj"] * Tai["Ak"] * Vijab["JoAh"] * Rabij["feJi"];
  ( + 0.5  ) * Tai["do"] * Tai["hk"] * Tai["Aj"] * Vijab["JoAh"] * Rabij["feJi"];

  ( - 0.5  ) * Tai["fo"] * Tai["ep"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bdjk"];
  ( + 0.5  ) * Tai["eo"] * Tai["fp"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bdjk"];
  ( + 0.5  ) * Tai["fo"] * Tai["dp"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bejk"];
  ( - 0.5  ) * Tai["do"] * Tai["fp"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bejk"];
  ( - 0.5  ) * Tai["eo"] * Tai["dp"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bfjk"];
  ( + 0.5  ) * Tai["do"] * Tai["ep"] * Tai["Ai"] * Vijab["poBA"] * Rabij["Bfjk"];
  ( + 0.5  ) * Tai["fo"] * Tai["ep"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Bdik"];
  ( - 0.5  ) * Tai["eo"] * Tai["fp"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Bdik"];
  ( - 0.5  ) * Tai["fo"] * Tai["dp"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Beik"];
  ( + 0.5  ) * Tai["do"] * Tai["fp"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Beik"];
  ( + 0.5  ) * Tai["eo"] * Tai["dp"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Bfik"];
  ( - 0.5  ) * Tai["do"] * Tai["ep"] * Tai["Aj"] * Vijab["poBA"] * Rabij["Bfik"];
  ( + 0.5  ) * Tai["fo"] * Tai["ep"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Bdji"];
  ( - 0.5  ) * Tai["eo"] * Tai["fp"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Bdji"];
  ( - 0.5  ) * Tai["fo"] * Tai["dp"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Beji"];
  ( + 0.5  ) * Tai["do"] * Tai["fp"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Beji"];
  ( + 0.5  ) * Tai["eo"] * Tai["dp"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Bfji"];
  ( - 0.5  ) * Tai["do"] * Tai["ep"] * Tai["Ak"] * Vijab["poBA"] * Rabij["Bfji"];

  ( - 0.5  ) * Tai["gi"] * Tai["hj"] * Tabij["efIk"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 0.5  ) * Tai["gj"] * Tai["hi"] * Tabij["efIk"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 0.5  ) * Tai["gi"] * Tai["hj"] * Tabij["fdIk"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 0.5  ) * Tai["gj"] * Tai["hi"] * Tabij["fdIk"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 0.5  ) * Tai["gi"] * Tai["hj"] * Tabij["deIk"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 0.5  ) * Tai["gj"] * Tai["hi"] * Tabij["deIk"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 0.5  ) * Tai["gi"] * Tai["hk"] * Tabij["efIj"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 0.5  ) * Tai["gk"] * Tai["hi"] * Tabij["efIj"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 0.5  ) * Tai["gi"] * Tai["hk"] * Tabij["fdIj"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 0.5  ) * Tai["gk"] * Tai["hi"] * Tabij["fdIj"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 0.5  ) * Tai["gi"] * Tai["hk"] * Tabij["deIj"] * Vijab["JIhg"] * Rai["fJ"];
  ( - 0.5  ) * Tai["gk"] * Tai["hi"] * Tabij["deIj"] * Vijab["JIhg"] * Rai["fJ"];
  ( - 0.5  ) * Tai["gj"] * Tai["hk"] * Tabij["efIi"] * Vijab["JIhg"] * Rai["dJ"];
  ( + 0.5  ) * Tai["gk"] * Tai["hj"] * Tabij["efIi"] * Vijab["JIhg"] * Rai["dJ"];
  ( - 0.5  ) * Tai["gj"] * Tai["hk"] * Tabij["fdIi"] * Vijab["JIhg"] * Rai["eJ"];
  ( + 0.5  ) * Tai["gk"] * Tai["hj"] * Tabij["fdIi"] * Vijab["JIhg"] * Rai["eJ"];
  ( - 0.5  ) * Tai["gj"] * Tai["hk"] * Tabij["deIi"] * Vijab["JIhg"] * Rai["fJ"];
  ( + 0.5  ) * Tai["gk"] * Tai["hj"] * Tabij["deIi"] * Vijab["JIhg"] * Rai["fJ"];

  ( + 1.0  ) * Tai["fo"] * Tai["hi"] * Tabij["Aejk"] * Vijab["JoAh"] * Rai["dJ"];
  ( - 1.0  ) * Tai["fo"] * Tai["hi"] * Tabij["Adjk"] * Vijab["JoAh"] * Rai["eJ"];
  ( - 1.0  ) * Tai["eo"] * Tai["hi"] * Tabij["Afjk"] * Vijab["JoAh"] * Rai["dJ"];
  ( + 1.0  ) * Tai["do"] * Tai["hi"] * Tabij["Afjk"] * Vijab["JoAh"] * Rai["eJ"];
  ( + 1.0  ) * Tai["eo"] * Tai["hi"] * Tabij["Adjk"] * Vijab["JoAh"] * Rai["fJ"];
  ( - 1.0  ) * Tai["do"] * Tai["hi"] * Tabij["Aejk"] * Vijab["JoAh"] * Rai["fJ"];
  ( - 1.0  ) * Tai["fo"] * Tai["hj"] * Tabij["Aeik"] * Vijab["JoAh"] * Rai["dJ"];
  ( + 1.0  ) * Tai["fo"] * Tai["hj"] * Tabij["Adik"] * Vijab["JoAh"] * Rai["eJ"];
  ( + 1.0  ) * Tai["eo"] * Tai["hj"] * Tabij["Afik"] * Vijab["JoAh"] * Rai["dJ"];
  ( - 1.0  ) * Tai["do"] * Tai["hj"] * Tabij["Afik"] * Vijab["JoAh"] * Rai["eJ"];
  ( - 1.0  ) * Tai["eo"] * Tai["hj"] * Tabij["Adik"] * Vijab["JoAh"] * Rai["fJ"];
  ( + 1.0  ) * Tai["do"] * Tai["hj"] * Tabij["Aeik"] * Vijab["JoAh"] * Rai["fJ"];
  ( - 1.0  ) * Tai["fo"] * Tai["hk"] * Tabij["Aeji"] * Vijab["JoAh"] * Rai["dJ"];
  ( + 1.0  ) * Tai["fo"] * Tai["hk"] * Tabij["Adji"] * Vijab["JoAh"] * Rai["eJ"];
  ( + 1.0  ) * Tai["eo"] * Tai["hk"] * Tabij["Afji"] * Vijab["JoAh"] * Rai["dJ"];
  ( - 1.0  ) * Tai["do"] * Tai["hk"] * Tabij["Afji"] * Vijab["JoAh"] * Rai["eJ"];
  ( - 1.0  ) * Tai["eo"] * Tai["hk"] * Tabij["Adji"] * Vijab["JoAh"] * Rai["fJ"];
  ( + 1.0  ) * Tai["do"] * Tai["hk"] * Tabij["Aeji"] * Vijab["JoAh"] * Rai["fJ"];

  ( - 1.0  ) * Tai["fo"] * Tai["hi"] * Tabij["deIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( + 1.0  ) * Tai["eo"] * Tai["hi"] * Tabij["dfIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( + 1.0  ) * Tai["do"] * Tai["hi"] * Tabij["feIj"] * Vijab["IoBh"] * Rai["Bk"];
  ( + 1.0  ) * Tai["fo"] * Tai["hi"] * Tabij["deIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tai["eo"] * Tai["hi"] * Tabij["dfIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tai["do"] * Tai["hi"] * Tabij["feIk"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tai["fo"] * Tai["hj"] * Tabij["deIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tai["eo"] * Tai["hj"] * Tabij["dfIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tai["do"] * Tai["hj"] * Tabij["feIi"] * Vijab["IoBh"] * Rai["Bk"];
  ( - 1.0  ) * Tai["fo"] * Tai["hk"] * Tabij["deIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tai["eo"] * Tai["hk"] * Tabij["dfIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( + 1.0  ) * Tai["do"] * Tai["hk"] * Tabij["feIi"] * Vijab["IoBh"] * Rai["Bj"];
  ( - 1.0  ) * Tai["fo"] * Tai["hj"] * Tabij["deIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( + 1.0  ) * Tai["eo"] * Tai["hj"] * Tabij["dfIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( + 1.0  ) * Tai["do"] * Tai["hj"] * Tabij["feIk"] * Vijab["IoBh"] * Rai["Bi"];
  ( + 1.0  ) * Tai["fo"] * Tai["hk"] * Tabij["deIj"] * Vijab["IoBh"] * Rai["Bi"];
  ( - 1.0  ) * Tai["eo"] * Tai["hk"] * Tabij["dfIj"] * Vijab["IoBh"] * Rai["Bi"];
  ( - 1.0  ) * Tai["do"] * Tai["hk"] * Tabij["feIj"] * Vijab["IoBh"] * Rai["Bi"];

  ( + 0.5  ) * Tai["fo"] * Tai["ep"] * Tabij["Adij"] * Vijab["poBA"] * Rai["Bk"];
  ( - 0.5  ) * Tai["eo"] * Tai["fp"] * Tabij["Adij"] * Vijab["poBA"] * Rai["Bk"];
  ( - 0.5  ) * Tai["fo"] * Tai["dp"] * Tabij["Aeij"] * Vijab["poBA"] * Rai["Bk"];
  ( + 0.5  ) * Tai["do"] * Tai["fp"] * Tabij["Aeij"] * Vijab["poBA"] * Rai["Bk"];
  ( + 0.5  ) * Tai["eo"] * Tai["dp"] * Tabij["Afij"] * Vijab["poBA"] * Rai["Bk"];
  ( - 0.5  ) * Tai["do"] * Tai["ep"] * Tabij["Afij"] * Vijab["poBA"] * Rai["Bk"];
  ( - 0.5  ) * Tai["fo"] * Tai["ep"] * Tabij["Adik"] * Vijab["poBA"] * Rai["Bj"];
  ( + 0.5  ) * Tai["eo"] * Tai["fp"] * Tabij["Adik"] * Vijab["poBA"] * Rai["Bj"];
  ( + 0.5  ) * Tai["fo"] * Tai["dp"] * Tabij["Aeik"] * Vijab["poBA"] * Rai["Bj"];
  ( - 0.5  ) * Tai["do"] * Tai["fp"] * Tabij["Aeik"] * Vijab["poBA"] * Rai["Bj"];
  ( - 0.5  ) * Tai["eo"] * Tai["dp"] * Tabij["Afik"] * Vijab["poBA"] * Rai["Bj"];
  ( + 0.5  ) * Tai["do"] * Tai["ep"] * Tabij["Afik"] * Vijab["poBA"] * Rai["Bj"];
  ( - 0.5  ) * Tai["fo"] * Tai["ep"] * Tabij["Adkj"] * Vijab["poBA"] * Rai["Bi"];
  ( + 0.5  ) * Tai["eo"] * Tai["fp"] * Tabij["Adkj"] * Vijab["poBA"] * Rai["Bi"];
  ( + 0.5  ) * Tai["fo"] * Tai["dp"] * Tabij["Aekj"] * Vijab["poBA"] * Rai["Bi"];
  ( - 0.5  ) * Tai["do"] * Tai["fp"] * Tabij["Aekj"] * Vijab["poBA"] * Rai["Bi"];
  ( - 0.5  ) * Tai["eo"] * Tai["dp"] * Tabij["Afkj"] * Vijab["poBA"] * Rai["Bi"];
  ( + 0.5  ) * Tai["do"] * Tai["ep"] * Tabij["Afkj"] * Vijab["poBA"] * Rai["Bi"];

}

// instantiate
template
class CcsdSimilarityTransformedHamiltonian<complex>;
template
class CcsdSimilarityTransformedHamiltonian<double>;
