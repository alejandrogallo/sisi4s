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
#include <algorithms/SimilarityTransformedHamiltonian.hpp>
#include <algorithms/StantonIntermediatesUCCSD.hpp>
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

#ifdef DEBUG
#define ST_DEBUG(msg) \
  LOG(1, "debug:STHam:") << __LINE__ << ":" << \
  "\x1b[33m" << msg << "\x1b[0m" << std::endl;
#else
#define ST_DEBUG(msg)
#endif

using namespace cc4s;

template <typename F>
PTR(CTF::Tensor<F>) SimilarityTransformedHamiltonian<F>::getTauABIJ() {

  if (Tau_abij) {
    return Tau_abij;
  }

  LOG(1, getAbbreviation())
    << "Building Tau_abij from Tai and Tabij"
    << std::endl;

  Tau_abij = NEW(CTF::Tensor<F>, *Tabij);
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

  return Tau_abij;

}

template <typename F>
SDFockVector<F> SimilarityTransformedHamiltonian<F>::rightApply(
  SDFockVector<F> &R
) {
  return useRightApplyIntermediates ?
    rightApplyIntermediates(R) : rightApplyHirata(R);
}

template <typename F>
SDFockVector<F> SimilarityTransformedHamiltonian<F>::rightApply_CCSD_IP(
  SDFockVector<F> &R
) {
  return useRightApplyIntermediates ?
    rightApplyIntermediates_CCSD_IP(R) : rightApplyHirata_CCSD_IP(R);
}


template <typename F>
SDFockVector<F>
SimilarityTransformedHamiltonian<F>::rightApplyIntermediates_CCSD_IP(
  SDFockVector<F> &R
) {
  SDFockVector<F> HR(R);
  PTR(CTF::Tensor<F>) Ri(R.get(0) );
  PTR(CTF::Tensor<F>) Raij(R.get(1) );
  PTR(CTF::Tensor<F>) HRi(HR.get(0) );
  PTR(CTF::Tensor<F>) HRaij(HR.get(1) );

  // For singles
  Wij = getIJ();
  Wia = getIA();
  Wijka = getIJKA();

  // For doubles
  Wab = getAB();
  Wiajk = getIAJK();
  Wiabj = getIABJ();
  Wijkl = getIJKL();

  (*HRi)["i"]  = 0.0;
  (*HRi)["i"] += (-1.0) * (*Wij)["mi"] * (*Ri)["m"];
  (*HRi)["i"] += (-1.0) * (*Wia)["mc"] * (*Raij)["cim"];
  (*HRi)["i"] += (-0.5) * (*Wijka)["nmic"] * (*Raij)["cmn"];

  (*HRaij)["aij"]  = 0.0;
  (*HRaij)["aij"] += (-1.0) * (*Wiajk)["maji"] * (*Ri)["m"];
  (*HRaij)["aij"] += (+1.0) * (*Wab)["ac"] * (*Raij)["cij"];
  // we have to antisymmetrize here
  (*HRaij)["aij"] += (-1.0) * (*Wij)["mi"] * (*Raij)["amj"];
  (*HRaij)["aij"] += (+1.0) * (*Wij)["mj"] * (*Raij)["ami"];
  // also antisymmetrize
  (*HRaij)["aij"] += (+1.0) * (*Wiabj)["maci"] * (*Raij)["cmj"];
  (*HRaij)["aij"] += (-1.0) * (*Wiabj)["macj"] * (*Raij)["cmi"];

  (*HRaij)["aij"] += (+0.5) * (*Wijkl)["mnij"] * (*Raij)["amn"];

  // three body term
  // WHHPHPH = Wijakbl := Vijbc * Tcakl (page 354 Shavitt 10.88)
  // (*Wijakbl)["nmaicj"] * (*Raij)["cmn"];
  // (*Wijakbl)["nmaicj"] = (*Vijab)["mnce"] * (*Tabij)["eaij"];

  //(*HRaij)["aij"] += (-0.5) * (*Vijab)["mnce"] * (*Tabij)["eaij"] * (*Raij)["cmn"];
  // Change order for preformance
  // TODO: Really review this little puppy
  (*HRaij)["aij"] += (-0.5) * (*Tabij)["eaji"] * (*Vijab)["mnce"] * (*Raij)["cmn"];

  return HR;
}

template <typename F>
SDFockVector<F>
SimilarityTransformedHamiltonian<F>::rightApplyHirata_CCSD_IP(
  SDFockVector<F> &R
) {

  SDFockVector<F> HR(R);
  PTR(CTF::Tensor<F>) Ri(R.get(0) );
  PTR(CTF::Tensor<F>) Raij(R.get(1) );
  PTR(CTF::Tensor<F>) HRi(HR.get(0) );
  PTR(CTF::Tensor<F>) HRaij(HR.get(1) );

  (*HRi)["i"] = 0.0;
  if (Fia) {
    (*HRi)["i"] += (+ 1.0) * (*Fia)["kc"] * (*Raij)["cki"];
    (*HRi)["i"] += (- 1.0) * (*Fia)["kc"] * (*Tai)["ci"] * (*Ri)["k"];
  }
  (*HRi)["i"] += (- 1.0) * (*Fij)["ki"] * (*Ri)["k"];
  (*HRi)["i"] += (+ 0.5) * (*Vijka)["klid"] * (*Raij)["dkl"];
  (*HRi)["i"] += (- 1.0) * (*Tai)["bl"] * (*Vijka)["mlib"] * (*Ri)["m"];
  (*HRi)["i"] += (- 0.5) * (*Tai)["bi"] * (*Vijab)["lmeb"] * (*Raij)["elm"];
  (*HRi)["i"] += (+ 1.0) * (*Tai)["bl"] * (*Vijab)["mleb"] * (*Raij)["emi"];
  (*HRi)["i"] += (+ 0.5) * (*Tabij)["bcmi"] * (*Vijab)["nmbc"] * (*Ri)["n"];
  (*HRi)["i"] += (+ 1.0) * (*Tai)["bi"] * (*Tai)["cm"] * (*Vijab)["nmcb"] * (*Ri)["n"];

  (*HRaij)["bij"] = 0.0;
  if (Fia) {
    //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bij=>bji"]) * (*Fia)["me"] * (*Tai)["ei"] * (*Raij)["bmj"];
    (*HRaij)["bij"] += (- 1.0) * (*Fia)["me"] * (*Tai)["ei"] * (*Raij)["bmj"];
    (*HRaij)["bij"] += (+ 1.0) * (*Fia)["me"] * (*Tai)["ej"] * (*Raij)["bmi"];
    (*HRaij)["bij"] += (- 1.0) * (*Fia)["me"] * (*Tai)["bm"] * (*Raij)["eij"];
    (*HRaij)["bij"] += (+ 1.0) * (*Fia)["me"] * (*Tabij)["ebij"] * (*Ri)["m"];
  }
  (*HRaij)["bij"] += (+ 1.0) * (*Viajk)["mbij"] * (*Ri)["m"];

  //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bij=>bji"]) * (*Fij)["mi"] * (*Raij)["bmj"];
  (*HRaij)["bij"] += (- 1.0) * (*Fij)["mi"] * (*Raij)["bmj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Fij)["mj"] * (*Raij)["bmi"];

  (*HRaij)["bij"] += (+ 1.0) * (*Fab)["bd"] * (*Raij)["dij"];
  (*HRaij)["bij"] += (+ 0.5) * (*Vijkl)["mnij"] * (*Raij)["bmn"];

  //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bij=>bji"]) * (*Viajb)["mbie"] * (*Raij)["emj"];
  (*HRaij)["bij"] += (- 1.0) * (*Viajb)["mbie"] * (*Raij)["emj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Viajb)["mbje"] * (*Raij)["emi"];

  (*HRaij)["bij"] += (- 1.0) * (*Tai)["bm"] * (*Vijkl)["nmij"] * (*Ri)["n"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bji=>bij"]) * (*Tai)["dj"] * (*Viajb)["nbid"] * (*Ri)["n"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["dj"] * (*Viajb)["nbid"] * (*Ri)["n"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["di"] * (*Viajb)["nbjd"] * (*Ri)["n"];

  //(*HRaij)["bij"] += (+ 0.5 - 0.5 * P["bji=>bij"]) * (*Tai)["dj"] * (*Vijka)["noid"] * (*Raij)["bno"];
  (*HRaij)["bij"] += (+ 0.5) * (*Tai)["dj"] * (*Vijka)["noid"] * (*Raij)["bno"];
  (*HRaij)["bij"] += (- 0.5) * (*Tai)["di"] * (*Vijka)["nojd"] * (*Raij)["bno"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bij=>bji"]) * (*Tai)["bm"] * (*Vijka)["nmif"] * (*Raij)["fnj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["bm"] * (*Vijka)["nmif"] * (*Raij)["fnj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["bm"] * (*Vijka)["nmjf"] * (*Raij)["fni"];

  //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bij=>bji"]) * (*Tai)["dn"] * (*Vijka)["onid"] * (*Raij)["boj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["dn"] * (*Vijka)["onid"] * (*Raij)["boj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["dn"] * (*Vijka)["onjd"] * (*Raij)["boi"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bij=>bji"]) * (*Tai)["di"] * (*Viabc)["nbfd"] * (*Raij)["fnj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["di"] * (*Viabc)["nbfd"] * (*Raij)["fnj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["dj"] * (*Viabc)["nbfd"] * (*Raij)["fni"];

  (*HRaij)["bij"] += (- 1.0) * (*Tai)["dn"] * (*Viabc)["nbfd"] * (*Raij)["fij"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bji=>bij"]) * (*Tabij)["dbnj"] * (*Vijka)["onid"] * (*Ri)["o"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tabij)["dbnj"] * (*Vijka)["onid"] * (*Ri)["o"];
  (*HRaij)["bij"] += (- 1.0) * (*Tabij)["dbni"] * (*Vijka)["onjd"] * (*Ri)["o"];


  (*HRaij)["bij"] += (+ 0.5) * (*Tabij)["deij"] * (*Viabc)["obde"] * (*Ri)["o"];
  (*HRaij)["bij"] += (+ 0.5) * (*Tabij)["dbij"] * (*Vijab)["nogd"] * (*Raij)["gno"];
  (*HRaij)["bij"] += (+ 0.25) * (*Tabij)["deij"] * (*Vijab)["oIde"] * (*Raij)["boI"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bij=>bji"]) * (*Tabij)["dbni"] * (*Vijab)["ongd"] * (*Raij)["goj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tabij)["dbni"] * (*Vijab)["ongd"] * (*Raij)["goj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tabij)["dbnj"] * (*Vijab)["ongd"] * (*Raij)["goi"];

  //(*HRaij)["bij"] += (+ 0.5 - 0.5 * P["bij=>bji"]) * (*Tabij)["deoi"] * (*Vijab)["Iode"] * (*Raij)["bIj"];
  (*HRaij)["bij"] += (+ 0.5) * (*Tabij)["deoi"] * (*Vijab)["Iode"] * (*Raij)["bIj"];
  (*HRaij)["bij"] += (- 0.5) * (*Tabij)["deoj"] * (*Vijab)["Iode"] * (*Raij)["bIi"];


  (*HRaij)["bij"] += (+ 0.5) * (*Tabij)["dbno"] * (*Vijab)["nogd"] * (*Raij)["gij"];

  //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bji=>bij"]) * (*Tai)["bm"] * (*Tai)["ej"] * (*Vijka)["omie"] * (*Ri)["o"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["bm"] * (*Tai)["ej"] * (*Vijka)["omie"] * (*Ri)["o"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["bm"] * (*Tai)["ei"] * (*Vijka)["omje"] * (*Ri)["o"];

  //(*HRaij)["bij"] += (- 0.5 + 0.5 * P["bij=>bji"]) * (*Tai)["di"] * (*Tai)["ej"] * (*Viabc)["obed"] * (*Ri)["o"];
  (*HRaij)["bij"] += (- 0.5) * (*Tai)["di"] * (*Tai)["ej"] * (*Viabc)["obed"] * (*Ri)["o"];
  (*HRaij)["bij"] += (+ 0.5) * (*Tai)["dj"] * (*Tai)["ei"] * (*Viabc)["obed"] * (*Ri)["o"];

  //(*HRaij)["bij"] += (- 0.25 + 0.25 * P["bij=>bji"]) * (*Tai)["di"] * (*Tai)["ej"] * (*Vijab)["oIed"] * (*Raij)["boI"];
  (*HRaij)["bij"] += (- 0.25) * (*Tai)["di"] * (*Tai)["ej"] * (*Vijab)["oIed"] * (*Raij)["boI"];
  (*HRaij)["bij"] += (+ 0.25) * (*Tai)["dj"] * (*Tai)["ei"] * (*Vijab)["oIed"] * (*Raij)["boI"];

  //(*HRaij)["bij"] += (- 1.0 + 1.0 * P["bij=>bji"]) * (*Tai)["bm"] * (*Tai)["ei"] * (*Vijab)["omge"] * (*Raij)["goj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["bm"] * (*Tai)["ei"] * (*Vijab)["omge"] * (*Raij)["goj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["bm"] * (*Tai)["ej"] * (*Vijab)["omge"] * (*Raij)["goi"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bij=>bji"]) * (*Tai)["di"] * (*Tai)["eo"] * (*Vijab)["Ioed"] * (*Raij)["bIj"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["di"] * (*Tai)["eo"] * (*Vijab)["Ioed"] * (*Raij)["bIj"];
  (*HRaij)["bij"] += (- 1.0) * (*Tai)["dj"] * (*Tai)["eo"] * (*Vijab)["Ioed"] * (*Raij)["bIi"];


  (*HRaij)["bij"] += (+ 1.0) * (*Tai)["bm"] * (*Tai)["eo"] * (*Vijab)["omge"] * (*Raij)["gij"];

  //(*HRaij)["bij"] += (+ 1.0 - 1.0 * P["bji=>bij"]) * (*Tabij)["dbnj"] * (*Tai)["fi"] * (*Vijab)["Infd"] * (*Ri)["I"];
  (*HRaij)["bij"] += (+ 1.0) * (*Tabij)["dbnj"] * (*Tai)["fi"] * (*Vijab)["Infd"] * (*Ri)["I"];
  (*HRaij)["bij"] += (- 1.0) * (*Tabij)["dbni"] * (*Tai)["fj"] * (*Vijab)["Infd"] * (*Ri)["I"];


  (*HRaij)["bij"] += (- 0.5) * (*Tabij)["deij"] * (*Tai)["bo"] * (*Vijab)["Iode"] * (*Ri)["I"];
  (*HRaij)["bij"] += (- 1.0) * (*Tabij)["dbij"] * (*Tai)["eo"] * (*Vijab)["Ioed"] * (*Ri)["I"];

  //(*HRaij)["bij"] += (+ 0.5 - 0.5 * P["bij=>bji"]) * (*Tai)["bm"] * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["Imfe"] * (*Ri)["I"];
  (*HRaij)["bij"] += (+ 0.5) * (*Tai)["bm"] * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["Imfe"] * (*Ri)["I"];
  (*HRaij)["bij"] += (- 0.5) * (*Tai)["bm"] * (*Tai)["ej"] * (*Tai)["fi"] * (*Vijab)["Imfe"] * (*Ri)["I"];


  return HR;;

}


template <typename F>
SDFockVector<F> SimilarityTransformedHamiltonian<F>::rightApply_CCSD_EA(
  SDFockVector<F> &R
) {

  SDFockVector<F> HR(R);
  PTR(CTF::Tensor<F>) Ra(R.get(0) );
  PTR(CTF::Tensor<F>) Rabi(R.get(1) );
  PTR(CTF::Tensor<F>) HRa(HR.get(0) );
  PTR(CTF::Tensor<F>) HRabi(HR.get(1) );

  (*HRa)["a"]  = 0.0;
  (*HRa)["a"] += (+ 1.0) * (*Fab)["ab"] * (*Ra)["b"];
  (*HRa)["a"] += (+ 1.0) * (*Fia)["kc"] * (*Rabi)["cak"];
  (*HRa)["a"] += (+ 0.5) * (*Viabc)["kacd"] * (*Rabi)["cdk"];
  (*HRa)["a"] += (- 1.0) * (*Fia)["kc"] * (*Tai)["ak"] * (*Ra)["c"];
  (*HRa)["a"] += (- 1.0) * (*Tai)["bl"] * (*Viabc)["ladb"] * (*Ra)["d"];
  (*HRa)["a"] += (- 0.5) * (*Tai)["ak"] * (*Vijab)["lkde"] * (*Rabi)["del"];
  (*HRa)["a"] += (+ 1.0) * (*Tai)["bl"] * (*Vijab)["mleb"] * (*Rabi)["eam"];
  (*HRa)["a"] += (+ 0.5) * (*Tabij)["balm"] * (*Vijab)["lmeb"] * (*Ra)["e"];
  (*HRa)["a"] += (+ 1.0) * (*Tai)["ak"] * (*Tai)["cm"] * (*Vijab)["mkec"] * (*Ra)["e"];

  (*HRabi)["bci"]  = 0.0;
  (*HRabi)["bci"] += (+ 1.0) * (*Vabic)["bcid"] * (*Ra)["d"];
  (*HRabi)["bci"] += (- 1.0) * (*Fij)["mi"] * (*Rabi)["bcm"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["cbi=>bci"]) * (*Fab)["cd"] * (*Rabi)["dbi"];
  (*HRabi)["bci"] += (- 1.0) * (*Fab)["cd"] * (*Rabi)["dbi"];
  (*HRabi)["bci"] += (+ 1.0) * (*Fab)["bd"] * (*Rabi)["dci"];

  //(*HRabi)["bci"] += (+ 1.0 - 1.0 * P["cbi=>bci"]) * (*Viajb)["mcie"] * (*Rabi)["ebm"];
  (*HRabi)["bci"] += (+ 1.0) * (*Viajb)["mcie"] * (*Rabi)["ebm"];
  (*HRabi)["bci"] += (- 1.0) * (*Viajb)["mbie"] * (*Rabi)["ecm"];

  (*HRabi)["bci"] += (+ 0.5) * (*Vabcd)["bcde"] * (*Rabi)["dei"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["bci=>cbi"]) * (*Tai)["bm"] * (*Viajb)["mcie"] * (*Ra)["e"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["bm"] * (*Viajb)["mcie"] * (*Ra)["e"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["cm"] * (*Viajb)["mbie"] * (*Ra)["e"];

  (*HRabi)["bci"] += (- 1.0) * (*Tai)["di"] * (*Vabcd)["bced"] * (*Ra)["e"];
  (*HRabi)["bci"] += (- 1.0) * (*Fia)["me"] * (*Tai)["ei"] * (*Rabi)["bcm"];

  //(*HRabi)["bci"] += (+ 1.0 - 1.0 * P["cbi=>bci"]) * (*Fia)["me"] * (*Tai)["cm"] * (*Rabi)["ebi"];
  (*HRabi)["bci"] += (+ 1.0) * (*Fia)["me"] * (*Tai)["cm"] * (*Rabi)["ebi"];
  (*HRabi)["bci"] += (- 1.0) * (*Fia)["me"] * (*Tai)["bm"] * (*Rabi)["eci"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Vijka)["nmif"] * (*Rabi)["fbn"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["cm"] * (*Vijka)["nmif"] * (*Rabi)["fbn"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["bm"] * (*Vijka)["nmif"] * (*Rabi)["fcn"];

  (*HRabi)["bci"] += (- 1.0) * (*Tai)["dn"] * (*Vijka)["onid"] * (*Rabi)["bco"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["cbi=>bci"]) * (*Tai)["di"] * (*Viabc)["ncfd"] * (*Rabi)["fbn"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["di"] * (*Viabc)["ncfd"] * (*Rabi)["fbn"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["di"] * (*Viabc)["nbfd"] * (*Rabi)["fcn"];

  //(*HRabi)["bci"] += (- 0.5 + 0.5 * P["bci=>cbi"]) * (*Tai)["bm"] * (*Viabc)["mcef"] * (*Rabi)["efi"];
  (*HRabi)["bci"] += (- 0.5) * (*Tai)["bm"] * (*Viabc)["mcef"] * (*Rabi)["efi"];
  (*HRabi)["bci"] += (+ 0.5) * (*Tai)["cm"] * (*Viabc)["mbef"] * (*Rabi)["efi"];

  //(*HRabi)["bci"] += (+ 1.0 - 1.0 * P["cbi=>bci"]) * (*Tai)["dn"] * (*Viabc)["ncfd"] * (*Rabi)["fbi"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["dn"] * (*Viabc)["ncfd"] * (*Rabi)["fbi"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["dn"] * (*Viabc)["nbfd"] * (*Rabi)["fci"];

  (*HRabi)["bci"] += (+ 1.0) * (*Fia)["me"] * (*Tabij)["bcmi"] * (*Ra)["e"];
  (*HRabi)["bci"] += (+ 0.5) * (*Tabij)["bcmn"] * (*Vijka)["mnif"] * (*Ra)["f"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["bci=>cbi"]) * (*Tabij)["dbni"] * (*Viabc)["ncfd"] * (*Ra)["f"];
  (*HRabi)["bci"] += (- 1.0) * (*Tabij)["dbni"] * (*Viabc)["ncfd"] * (*Ra)["f"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tabij)["dcni"] * (*Viabc)["nbfd"] * (*Ra)["f"];

  (*HRabi)["bci"] += (+ 0.5) * (*Tabij)["bcmi"] * (*Vijab)["nmfg"] * (*Rabi)["fgn"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["cbi=>bci"]) * (*Tabij)["dcni"] * (*Vijab)["ongd"] * (*Rabi)["gbo"];
  (*HRabi)["bci"] += (- 1.0) * (*Tabij)["dcni"] * (*Vijab)["ongd"] * (*Rabi)["gbo"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tabij)["dbni"] * (*Vijab)["ongd"] * (*Rabi)["gco"];

  (*HRabi)["bci"] += (+ 0.5) * (*Tabij)["deoi"] * (*Vijab)["Iode"] * (*Rabi)["bcI"];
  (*HRabi)["bci"] += (+ 0.25) * (*Tabij)["bcmn"] * (*Vijab)["mnfg"] * (*Rabi)["fgi"];

  //(*HRabi)["bci"] += (- 0.5 + 0.5 * P["cbi=>bci"]) * (*Tabij)["dcno"] * (*Vijab)["nogd"] * (*Rabi)["gbi"];
  (*HRabi)["bci"] += (- 0.5) * (*Tabij)["dcno"] * (*Vijab)["nogd"] * (*Rabi)["gbi"];
  (*HRabi)["bci"] += (+ 0.5) * (*Tabij)["dbno"] * (*Vijab)["nogd"] * (*Rabi)["gci"];

  //(*HRabi)["bci"] += (+ 0.5 - 0.5 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Tai)["bn"] * (*Vijka)["nmif"] * (*Ra)["f"];
  (*HRabi)["bci"] += (+ 0.5) * (*Tai)["cm"] * (*Tai)["bn"] * (*Vijka)["nmif"] * (*Ra)["f"];
  (*HRabi)["bci"] += (- 0.5) * (*Tai)["bm"] * (*Tai)["cn"] * (*Vijka)["nmif"] * (*Ra)["f"];

  //(*HRabi)["bci"] += (+ 1.0 - 1.0 * P["bci=>cbi"]) * (*Tai)["bm"] * (*Tai)["ei"] * (*Viabc)["mcfe"] * (*Ra)["f"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["bm"] * (*Tai)["ei"] * (*Viabc)["mcfe"] * (*Ra)["f"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["cm"] * (*Tai)["ei"] * (*Viabc)["mbfe"] * (*Ra)["f"];

  //(*HRabi)["bci"] += (+ 1.0 - 1.0 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Tai)["ei"] * (*Vijab)["omge"] * (*Rabi)["gbo"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["cm"] * (*Tai)["ei"] * (*Vijab)["omge"] * (*Rabi)["gbo"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["bm"] * (*Tai)["ei"] * (*Vijab)["omge"] * (*Rabi)["gco"];

  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["di"] * (*Tai)["eo"] * (*Vijab)["Ioed"] * (*Rabi)["bcI"];

  //(*HRabi)["bci"] += (+ 0.25 - 0.25 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Tai)["bn"] * (*Vijab)["nmfg"] * (*Rabi)["fgi"];
  (*HRabi)["bci"] += (+ 0.25) * (*Tai)["cm"] * (*Tai)["bn"] * (*Vijab)["nmfg"] * (*Rabi)["fgi"];
  (*HRabi)["bci"] += (- 0.25) * (*Tai)["bm"] * (*Tai)["cn"] * (*Vijab)["nmfg"] * (*Rabi)["fgi"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Tai)["eo"] * (*Vijab)["omge"] * (*Rabi)["gbi"];
  (*HRabi)["bci"] += (- 1.0) * (*Tai)["cm"] * (*Tai)["eo"] * (*Vijab)["omge"] * (*Rabi)["gbi"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tai)["bm"] * (*Tai)["eo"] * (*Vijab)["omge"] * (*Rabi)["gci"];

  (*HRabi)["bci"] += (- 0.5) * (*Tabij)["bcmn"] * (*Tai)["fi"] * (*Vijab)["mngf"] * (*Ra)["g"];

  //(*HRabi)["bci"] += (- 1.0 + 1.0 * P["bci=>cbi"]) * (*Tabij)["dbni"] * (*Tai)["co"] * (*Vijab)["ongd"] * (*Ra)["g"];
  (*HRabi)["bci"] += (- 1.0) * (*Tabij)["dbni"] * (*Tai)["co"] * (*Vijab)["ongd"] * (*Ra)["g"];
  (*HRabi)["bci"] += (+ 1.0) * (*Tabij)["dcni"] * (*Tai)["bo"] * (*Vijab)["ongd"] * (*Ra)["g"];

  (*HRabi)["bci"] += (- 1.0) * (*Tabij)["bcmi"] * (*Tai)["eo"] * (*Vijab)["omge"] * (*Ra)["g"];

  //(*HRabi)["bci"] += (- 0.5 + 0.5 * P["cbi=>bci"]) * (*Tai)["cm"] * (*Tai)["bn"] * (*Tai)["fi"] * (*Vijab)["nmgf"] * (*Ra)["g"];
  (*HRabi)["bci"] += (- 0.5) * (*Tai)["cm"] * (*Tai)["bn"] * (*Tai)["fi"] * (*Vijab)["nmgf"] * (*Ra)["g"];
  (*HRabi)["bci"] += (+ 0.5) * (*Tai)["bm"] * (*Tai)["cn"] * (*Tai)["fi"] * (*Vijab)["nmgf"] * (*Ra)["g"];

  return HR;;

}

template <typename F>
SFockVector<F> SimilarityTransformedHamiltonian<F>::rightApplyHirata_RPA(
  SFockVector<F> &R
) {

  // This is only using Viajb and Vijab
  SFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );

  int ovvo[] = {No, Nv, Nv, No};
  int oo[] = {No, No};
  int vv[] = {Nv, Nv};
  int lens2[] = {NS, NS};
  int lens[] = {NS, NS, NS, NS};

  if (!Wiabj) {
    Wiabj = NEW(CTF::Tensor<F>, 4, ovvo, lens, *Cc4s::world, "Wiabj");
    ST_DEBUG("build Wiabj")
    (*Wiabj)["iabj"] = (*Tabij)["camj"] * (*Vijab)["micb"];
  }

  if (!Wab) {
    Wab = NEW(CTF::Tensor<F>, 2, vv, lens2, *Cc4s::world, "Wab");
    ST_DEBUG("build Wab")
    (*Wab)["ab"]  = (*Fab)["ab"];
    (*Wab)["ab"] += (- 0.5) * (*Vijab)["mnbe"] * (*Tabij)["aemn"];
  }

  if (!Wij) {
    Wij = NEW(CTF::Tensor<F>, 2, oo, lens2, *Cc4s::world, "Wij");
    ST_DEBUG("build Wij")
    (*Wij)["ij"]  = (*Fij)["ij"];
    (*Wij)["ij"] += (  0.5) * (*Vijab)["imef"] * (*Tabij)["efjm"];
  }

  ST_DEBUG("rightApplyHirata_RPA")

  // Contruct HR (one body part)
  // TODO: why "bi" not "ai"?
  (*HRai)["bi"]  = 0.0;

  //ST_DEBUG("WIJ")
  // WIJ =====================================================================
  //(*HRai)["bi"] += ( - 1.0 ) * (*Fij)["ki"] * (*Rai)["bk"];
  //(*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  (*HRai)["bi"] += ( - 1.0 ) * (*Wij)["ni"] * (*Rai)["bn"];
  if (Tai) {
    (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  }

  //ST_DEBUG("WAB")
  // WAB =====================================================================
  //(*HRai)["bi"] += ( + 1.0 ) * (*Fab)["bc"] * (*Rai)["ci"];
  //(*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];
  (*HRai)["bi"] += (*Wab)["bc"] * (*Rai)["ci"];
  if (Tai) {
    (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijab)["kmdf"] * (*Rai)["fi"];
  }

  //ST_DEBUG("WIABJ")
  // WIABJ ===================================================================
  //(*HRai)["bi"] += ( + 1.0 )  * (*Rai)["fm"] * (*Tabij)["cbli"] * (*Vijab)["lmcf"];
  (*HRai)["bi"] += ( + 1.0 )  * (*Wiabj)["mbfi"] * (*Rai)["fm"];
  if (Tai) {
    ST_DEBUG("T1")
    (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
  }

  return HR;

}

template <typename F>
SDFockVector<F> SimilarityTransformedHamiltonian<F>::rightApplyHirata(
  SDFockVector<F> &R
) {
  SDFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  ST_DEBUG("rightApplyHirata Ccsd")

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
SDFockVector<F> SimilarityTransformedHamiltonian<F>::rightApplyIntermediates(
  SDFockVector<F> &R
) {
  SDFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  ST_DEBUG("rightApplyIntermediates Ccsd")

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

  (*HRai)["ai"]  = 0.0;
  (*HRai)["ai"] += (- 1.0) * (*Wij)["li"] * (*Rai)["al"];
  (*HRai)["ai"] += (*Wab)["ad"] * (*Rai)["di"];
  (*HRai)["ai"] += (*Wiabj)["ladi"] * (*Rai)["dl"];

  (*HRai)["ai"] += (*Wia)["ld"] * (*Rabij)["adil"];

  (*HRai)["ai"] += ( - 0.5 ) * (*Wijka)["lmid"] * (*Rabij)["adlm"];
  (*HRai)["ai"] += (   0.5 ) * (*Waibc)["alde"] * (*Rabij)["deil"];

  ST_DEBUG("singles done")

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

  ST_DEBUG("doubles done")

  return HR;
}

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getIJ() {
  if (Wij) return Wij;

  Wij = NEW(CTF::Tensor<F>, *Fij);

  if (true) {
    // This is the second row in diagram 10.55 in [1]
    LOG(1, getAbbreviation()) << "Building Wij from Wia" << std::endl;
    Wia = getIA();
    (*Wij)["ij"]  = (*Fij)["ij"];
    (*Wij)["ij"] += (*Vijka)["imje"] * (*Tai)["em"];
    (*Wij)["ij"] += (*Wia)["ie"] * (*Tai)["ej"];
    (*Wij)["ij"] += (+ 0.5) * (*Vijab)["imef"] * (*Tabij)["efjm"];
  } else {
    Tau_abij = getTauABIJ();
    // This is the first row in diagram 10.55 in [1]
    LOG(1, getAbbreviation()) << "Building Wij" << std::endl;
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
SimilarityTransformedHamiltonian<F>::getAB() {
  if (Wab) return Wab;

  Wab   = NEW(CTF::Tensor<F>, *Fab);

  if (true) {
    //diagram (10.54) second line in [1]
    LOG(1, getAbbreviation()) << "Building Wab from Wia" << std::endl;
    Wia = getIA();
    (*Wab)["ab"]  = (*Fab)["ab"];
    (*Wab)["ab"] += (*Viabc)["mafb"] * (*Tai)["fm"];
    (*Wab)["ab"] += ( -1.0) * (*Wia)["mb"] * (*Tai)["am"];
    (*Wab)["ab"] += (- 0.5) * (*Vijab)["mnbe"] * (*Tabij)["aemn"];
  } else {
    Tau_abij = getTauABIJ();
    //diagram (10.54) first line in [1]
    LOG(1, getAbbreviation()) << "Building Wab" << std::endl;
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
SimilarityTransformedHamiltonian<F>::getAI() {
  if (Wai) return Wai;
  LOG(1, getAbbreviation()) << "Building Wai" << std::endl;

  int syms[] = {NS, NS};
  int ov[] = {Nv, No};
  CTF::Tensor<F> InitFai(2, ov, syms, *Cc4s::world, "InitFai");

  Wai = NEW(CTF::Tensor<F>, InitFai);

  (*Wai)["bi"] = 0.0;
  if (dressing == Dressing(CCSD)) {
    LOG(1, getAbbreviation()) << "Wai = 0 since CCSD" << std::endl;
    return Wai;
  }

  if (useStantonIntermediatesUCCSD()) {

    auto intermediates = getStantonIntermediatesUCCSD();
    (*Wai)["bi"] = (*intermediates->getRai())["bi"];
    (*Wai)["bi"] += (+ 1.0) * (*Fab)["bc"] * (*Tai)["ci"];
    (*Wai)["bi"] += (- 1.0) * (*Fij)["ki"] * (*Tai)["bk"];

  } else {

    // Equations from Hirata
    if (Fia) {
      (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["ib"];
      (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["kd"] * (*Tabij)["dbki"];
      (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
    }
    (*Wai)["bi"] += (+ 1.0) * (*Fab)["bc"] * (*Tai)["ci"];
    (*Wai)["bi"] += (- 1.0) * (*Fij)["ki"] * (*Tai)["bk"];
    (*Wai)["bi"] += (- 1.0) * (*Tai)["cl"] * (*Viajb)["lbic"];
    (*Wai)["bi"] += (+ 0.5) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
    (*Wai)["bi"] += (+ 0.5) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
    (*Wai)["bi"] += (- 1.0) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
    (*Wai)["bi"] += (- 1.0) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
    (*Wai)["bi"] += (- 0.5) * (*Tai)["fi"] * (*Tabij)["cblm"] * (*Vijab)["lmcf"];
    (*Wai)["bi"] += (- 0.5) * (*Tai)["bn"] * (*Tabij)["cdmi"] * (*Vijab)["mncd"];
    (*Wai)["bi"] += (+ 1.0) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
    (*Wai)["bi"] += (- 1.0) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];

  }
  if (Tabcijk) { // Begin Wai with Triples Tabcijk
    (*Wai)["bi"] += ( + 0.25 ) * (*Tabcijk)["cdbmni"] * (*Vijab)["mncd"];
  } // end Wai with Triples Tabcijk


  return Wai;
}

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getAI_RPA() {
  if (Wai) return Wai;
  LOG(1, getAbbreviation())
    << "Building Wai only with Vijab" << std::endl;

  int syms[] = {NS, NS};
  int ov[] = {Nv, No};
  CTF::Tensor<F> InitFai(2, ov, syms, *Cc4s::world, "InitFai");

  Wai = NEW(CTF::Tensor<F>, InitFai);

  (*Wai)["bi"] = 0.0;
  if (dressing == Dressing(RPA)) {
    LOG(1, getAbbreviation()) << "Wai = 0 since CCSD" << std::endl;
    return Wai;
  }

  // Equations from Hirata
  if (Fia) {
    (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["ib"];
    (*Wai)["bi"] += ( + 1.0  ) * (*Fia)["kd"] * (*Tabij)["dbki"];
    (*Wai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
  }
  (*Wai)["bi"] += (+ 1.0) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Wai)["bi"] += (- 1.0) * (*Fij)["ki"] * (*Tai)["bk"];

  (*Wai)["bi"] += (- 0.5) * (*Tai)["fi"] * (*Tabij)["cblm"] * (*Vijab)["lmcf"];
  (*Wai)["bi"] += (- 0.5) * (*Tai)["bn"] * (*Tabij)["cdmi"] * (*Vijab)["mncd"];
  (*Wai)["bi"] += (+ 1.0) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
  (*Wai)["bi"] += (- 1.0) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];

  if (Tabcijk) { // Begin Wai with Triples Tabcijk
    (*Wai)["bi"] += ( + 0.25 ) * (*Tabcijk)["cdbmni"] * (*Vijab)["mncd"];
  } // end Wai with Triples Tabcijk

  return Wai;
}


template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getIA() {
  if (Wia) return Wia;
  LOG(1, getAbbreviation()) << "Building Wia" << std::endl;

  int syms[] = {NS, NS};
  int ov[] = {No, Nv};
  CTF::Tensor<F> InitFia(2, ov, syms, *Cc4s::world, "InitFia");

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
SimilarityTransformedHamiltonian<F>::getIJAB() {
  if (Wijab) return Wijab;

  LOG(1, getAbbreviation()) << "Building Wijab = Vijab" << std::endl;
  Wijab = NEW(CTF::Tensor<F>, *Vijab);

  return Wijab;
}

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getABCD() {
  if (Wabcd) return Wabcd;
  LOG(1, getAbbreviation()) << "Building Wabcd" << std::endl;

  Tau_abij = getTauABIJ();
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
SimilarityTransformedHamiltonian<F>::getABCI() {
  if (Wabci) return Wabci;

  Wabci = NEW(CTF::Tensor<F>, *Vabci);

  bool wabciIntermediates(false);
  if (wabciIntermediates) {
    LOG(1, getAbbreviation()) << "Building Wabci from Wabcd and Wia"
                                         << std::endl;
    Tau_abij = getTauABIJ();
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
    LOG(1, getAbbreviation()) << "Building Wabci" << std::endl;
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
SimilarityTransformedHamiltonian<F>::getAIBC() {
  if (Waibc) return Waibc;
  LOG(1, getAbbreviation()) << "Building Waibc" << std::endl;

  Waibc = NEW(CTF::Tensor<F>, *Vaibc);

  (*Waibc)["aibc"]  = (*Vaibc)["aibc"];
  (*Waibc)["aibc"] += ( -1.0) * (*Vijab)["mibc"] * (*Tai)["am"];

  return Waibc;
}

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getIABJ() {
  if (Wiabj) return Wiabj;
  LOG(1, getAbbreviation()) << "Building Wiabj = Waijb" << std::endl;

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
SimilarityTransformedHamiltonian<F>::getIAJK() {
  if (Wiajk) return Wiajk;

  Wiajk = NEW(CTF::Tensor<F>, *Viajk);
  Wia = getIA();
  Wijkl = getIJKL();
  Tau_abij = getTauABIJ();

  LOG(1, getAbbreviation()) << "Building Wiajk from Wia and Wijkl" << std::endl;
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
SimilarityTransformedHamiltonian<F>::getIJKA() {
  if (Wijka) return Wijka;
  LOG(1, getAbbreviation()) << "Building Wijka" << std::endl;

  Wijka = NEW(CTF::Tensor<F>, *Vijka);

  //Taken directly from[2]
  (*Wijka)["jkia"]  = (*Vijka)["jkia"];
  (*Wijka)["jkia"] += (*Tai)["ei"] * (*Vijab)["jkea"];

  return Wijka;
}

template <typename F>
PTR(CTF::Tensor<F>)
SimilarityTransformedHamiltonian<F>::getIJKL() {
  if (Wijkl) return Wijkl;
  LOG(1, getAbbreviation()) << "Building Wijkl" << std::endl;

  Wijkl = NEW(CTF::Tensor<F>, *Vijkl);
  Tau_abij = getTauABIJ();

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
SDFockVector<F> SimilarityTransformedHamiltonian<F>::leftApplyIntermediates(
  SDFockVector<F> &L
){
  /*
  Equations from:

  @book{book:293288,
    author = {Isaiah Shavitt, Rodney J. Bartlett},
    edition = {1},
    isbn = {052181832X,9780521818322},
    publisher = {Cambridge University Press},
    series = {Cambridge Molecular Science},
    title = {Many-Body Methods in Chemistry and Physics: MBPT and Coupled-Cluster Theory},
    type = {book},
    year = {2009},
  }
  */
  SDFockVector<F> LH(L);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Lia( L.get(0) );
  PTR(CTF::Tensor<F>) Lijab( L.get(1) );
  PTR(CTF::Tensor<F>) LHia( LH.get(0) );
  PTR(CTF::Tensor<F>) LHijab( LH.get(1) );

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

  // Page 376
  (*LHia)["ia"]  = 0.0;
  // Diagram 1
  (*LHia)["ia"] += (*Wia)["ia"];

  // Diagram 3
  (*LHia)["ia"] += ( + 1.0  ) * (*Wab)["ea"] * (*Lia)["ie"];
  // Diagram 4
  (*LHia)["ia"] += ( - 1.0  ) * (*Wij)["ik"] * (*Lia)["ka"];
  // Diagram 5
  (*LHia)["ia"] += ( + 1.0  ) * (*Wiabj)["ieak"] * (*Lia)["ke"];
  // Diagram 18
  (*LHia)["ia"] += ( + 1.0  ) * (*Wabci)["edak"] * (*Lijab)["iked"];
  // Diagram 19
  (*LHia)["ia"] += ( - 1.0  ) * (*Wiajk)["iekl"] * (*Lijab)["klae"];
  // Diagram 26
  (*LHia)["ia"] +=
    ( - 1.0  ) * (*Waibc)["fida"] * (*Tabij)["edkl"] * (*Lijab)["klef"];
  // Diagram 27
  (*LHia)["ia"] +=
    ( - 1.0  ) * (*Wijka)["mila"] * (*Tabij)["edkm"] * (*Lijab)["kled"];

  // Page 377
  //Diagram 1
  (*LHijab)["ijab"] += ( + 1.0  ) * (*Vijab)["ijab"];
  //Diagram 2
  (*LHijab)["ijab"] += ( + 1.0  ) * (*LHia)["ia"] * (*Wia)["jb"];
  //Diagram 4
  (*LHijab)["ijab"] += ( + 1.0  ) * (*LHia)["ie"] * (*Waibc)["ejab"];
  // TODO: 5 8 9 10 11 12 23 24

  return LH;

}

template <typename F>
SDFockVector<F> SimilarityTransformedHamiltonian<F>::leftApplyHirata(
  SDFockVector<F> &L
) {
  SDFockVector<F> LH(L);
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
SDTFockVector<F> SimilarityTransformedHamiltonian<F>::rightApply(
  SDTFockVector<F> &R
) {
  return rightApplyHirata(R);
}

template <typename F>
SDTFockVector<F> SimilarityTransformedHamiltonian<F>::rightApplyHirata(
  SDTFockVector<F> &R
) {

  SDTFockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) Rabcijk( R.get(2) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );
  PTR(CTF::Tensor<F>) HRabcijk( HR.get(2) );

  { // keep CcsdR only in this scope
    SDFockVector<F> CcsdR(
        std::vector<PTR(CTF::Tensor<F>)>({Rai, Rabij}),
        std::vector<std::string>({"ai", "abij"})
    );

    SDFockVector<F> HCssdR(rightApply(CcsdR));

    (*HRai)["ai"] = (*HCssdR.get(0))["ai"];
    (*HRabij)["abij"] = (*HCssdR.get(1))["abij"];
  }

  ST_DEBUG("rightApplyHirata Ccsdt")

  //: BEGIN SINGLES
  (*HRai)["bi"] += ( + 0.25  ) * (*Vijab)["klef"] * (*Rabcijk)["feblki"];

  ST_DEBUG("singles done")

  //: BEGIN DOUBLES
  if (Fia) {
    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Rabcijk)["fcdmij"];
  }
  (*HRabij)["cdij"] += ( + 1.0  ) * (*Fab)["fg"] * (*Rabcijk)["gdeijk"];
  (*HRabij)["cdij"] += ( - 1.0  ) * (*Fab)["eg"] * (*Rabcijk)["gdfijk"];
  (*HRabij)["cdij"] += ( - 1.0  ) * (*Fab)["dg"] * (*Rabcijk)["gfeijk"];
  (*HRabij)["cdij"] += ( - 1.0  ) * (*Fij)["oi"] * (*Rabcijk)["defojk"];
  (*HRabij)["cdij"] += ( + 1.0  ) * (*Fij)["oj"] * (*Rabcijk)["defoik"];
  (*HRabij)["cdij"] += ( + 1.0  ) * (*Fij)["ok"] * (*Rabcijk)["defoji"];

  (*HRabij)["cdij"] += ( - 0.5  ) * (*Vijka)["mnig"] * (*Rabcijk)["gcdnmj"];
  (*HRabij)["cdij"] += ( + 0.5  ) * (*Vijka)["mnjg"] * (*Rabcijk)["gcdnmi"];

  (*HRabij)["cdij"] += ( + 0.5  ) * (*Viabc)["mdfg"] * (*Rabcijk)["gfcmij"];
  (*HRabij)["cdij"] += ( - 0.5  ) * (*Viabc)["mcfg"] * (*Rabcijk)["gfdmij"];

  if (dressing != NONE) {

    (*HRabij)["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Vijab)["nohe"] * (*Rabcijk)["hcdonj"];
    (*HRabij)["cdij"] += ( - 0.5  ) * (*Tai)["ej"] * (*Vijab)["nohe"] * (*Rabcijk)["hcdoni"];

    (*HRabij)["cdij"] += ( - 0.5  ) * (*Tai)["dm"] * (*Vijab)["nmgh"] * (*Rabcijk)["hgcnij"];
    (*HRabij)["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Vijab)["nmgh"] * (*Rabcijk)["hgdnij"];

    (*HRabij)["cdij"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijab)["onhe"] * (*Rabcijk)["hcdoij"];

    if (dressing == CCSDT) {

      (*HRabij)["cdij"] += ( + 1.0  ) * (*Tabcijk)["ecdnij"] * (*Vijab)["onhe"] * (*Rai)["ho"];

      (*HRabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["efdoij"] * (*Vijab)["poef"] * (*Rai)["cp"];
      (*HRabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["efcoij"] * (*Vijab)["poef"] * (*Rai)["dp"];

      (*HRabij)["cdij"] += ( - 0.5  ) * (*Tabcijk)["ecdnoi"] * (*Vijab)["nohe"] * (*Rai)["hj"];
      (*HRabij)["cdij"] += ( + 0.5  ) * (*Tabcijk)["ecdnoj"] * (*Vijab)["nohe"] * (*Rai)["hi"];

    } // Dressing(CCSD)

  } // Dressing(NONE)

  ST_DEBUG("doubles done")

  //: BEGIN TRIPLES
  if (Fia && dressing != NONE) {

    if (dressing == CCSDT) {
      (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hefijk"] * (*Rai)["do"];
      (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hfdijk"] * (*Rai)["eo"];
      (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["hdeijk"] * (*Rai)["fo"];

      (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defoij"] * (*Rai)["hk"];
      (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defoik"] * (*Rai)["hj"];
      (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabcijk)["defokj"] * (*Rai)["hi"];
    }

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tai)["hi"] * (*Rabcijk)["defojk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tai)["hj"] * (*Rabcijk)["defoik"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tai)["hk"] * (*Rabcijk)["defoji"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tai)["fo"] * (*Rabcijk)["hdeijk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tai)["eo"] * (*Rabcijk)["hdfijk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tai)["do"] * (*Rabcijk)["hfeijk"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["hfij"] * (*Rabij)["deok"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["heij"] * (*Rabij)["dfok"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["hdij"] * (*Rabij)["feok"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["hfik"] * (*Rabij)["deoj"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["heik"] * (*Rabij)["dfoj"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["hdik"] * (*Rabij)["feoj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["hfkj"] * (*Rabij)["deoi"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["hekj"] * (*Rabij)["dfoi"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["hdkj"] * (*Rabij)["feoi"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["efoi"] * (*Rabij)["hdjk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["fdoi"] * (*Rabij)["hejk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Fia)["oh"] * (*Tabij)["deoi"] * (*Rabij)["hfjk"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["efoj"] * (*Rabij)["hdik"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["fdoj"] * (*Rabij)["heik"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["deoj"] * (*Rabij)["hfik"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["efok"] * (*Rabij)["hdji"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["fdok"] * (*Rabij)["heji"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Fia)["oh"] * (*Tabij)["deok"] * (*Rabij)["hfji"];
  }

  ST_DEBUG("Vhphh * R2")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajk)["ofij"] * (*Rabij)["deok"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajk)["oeij"] * (*Rabij)["dfok"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajk)["odij"] * (*Rabij)["feok"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajk)["ofik"] * (*Rabij)["deoj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajk)["oeik"] * (*Rabij)["dfoj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajk)["odik"] * (*Rabij)["feoj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajk)["ofkj"] * (*Rabij)["deoi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajk)["oekj"] * (*Rabij)["dfoi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajk)["odkj"] * (*Rabij)["feoi"];

  ST_DEBUG("Vpphp * R2")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Vabic)["efig"] * (*Rabij)["gdjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Vabic)["fdig"] * (*Rabij)["gejk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Vabic)["deig"] * (*Rabij)["gfjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["efjg"] * (*Rabij)["gdik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["fdjg"] * (*Rabij)["geik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["dejg"] * (*Rabij)["gfik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["efkg"] * (*Rabij)["gdji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["fdkg"] * (*Rabij)["geji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Vabic)["dekg"] * (*Rabij)["gfji"];

  ST_DEBUG("Vhhhh * R3")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Vijkl)["opij"] * (*Rabcijk)["defpok"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Vijkl)["opik"] * (*Rabcijk)["defpoj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Vijkl)["opkj"] * (*Rabcijk)["defpoi"];

  ST_DEBUG("Vhphp * R3")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajb)["ofih"] * (*Rabcijk)["hdeojk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajb)["oeih"] * (*Rabcijk)["hdfojk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajb)["odih"] * (*Rabcijk)["hfeojk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajb)["ofjh"] * (*Rabcijk)["hdeoik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajb)["oejh"] * (*Rabcijk)["hdfoik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajb)["odjh"] * (*Rabcijk)["hfeoik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Viajb)["ofkh"] * (*Rabcijk)["hdeoji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajb)["oekh"] * (*Rabcijk)["hdfoji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Viajb)["odkh"] * (*Rabcijk)["hfeoji"];

  ST_DEBUG("Vpppp * R3")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Vabcd)["efgh"] * (*Rabcijk)["hgdijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Vabcd)["fdgh"] * (*Rabcijk)["hgeijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Vabcd)["degh"] * (*Rabcijk)["hgfijk"];

  if (dressing == NONE) {
    return HR;
  }

  ST_DEBUG("T1 * Vhhhh * R2")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Vijkl)["poij"] * (*Rabij)["depk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Vijkl)["poij"] * (*Rabij)["dfpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Vijkl)["poij"] * (*Rabij)["fepk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Vijkl)["poik"] * (*Rabij)["depj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Vijkl)["poik"] * (*Rabij)["dfpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Vijkl)["poik"] * (*Rabij)["fepj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Vijkl)["pokj"] * (*Rabij)["depi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Vijkl)["pokj"] * (*Rabij)["dfpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Vijkl)["pokj"] * (*Rabij)["fepi"];

  ST_DEBUG("T1 * Vhphp * R2")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Viajb)["pfig"] * (*Rabij)["depk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gj"] * (*Viajb)["peig"] * (*Rabij)["dfpk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gj"] * (*Viajb)["pdig"] * (*Rabij)["fepk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gk"] * (*Viajb)["pfig"] * (*Rabij)["depj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Viajb)["peig"] * (*Rabij)["dfpj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Viajb)["pdig"] * (*Rabij)["fepj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Viajb)["pfjg"] * (*Rabij)["depk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gi"] * (*Viajb)["pejg"] * (*Rabij)["dfpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gi"] * (*Viajb)["pdjg"] * (*Rabij)["fepk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gi"] * (*Viajb)["pfkg"] * (*Rabij)["depj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Viajb)["pekg"] * (*Rabij)["dfpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Viajb)["pdkg"] * (*Rabij)["fepj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Viajb)["pfjg"] * (*Rabij)["depi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gk"] * (*Viajb)["pejg"] * (*Rabij)["dfpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gk"] * (*Viajb)["pdjg"] * (*Rabij)["fepi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gj"] * (*Viajb)["pfkg"] * (*Rabij)["depi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Viajb)["pekg"] * (*Rabij)["dfpi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Viajb)["pdkg"] * (*Rabij)["fepi"];

  ST_DEBUG("T1 * Vhphp * R2")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Viajb)["ofih"] * (*Rabij)["hdjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Viajb)["ofih"] * (*Rabij)["hejk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Viajb)["oeih"] * (*Rabij)["hdjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Viajb)["odih"] * (*Rabij)["hejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Viajb)["oeih"] * (*Rabij)["hfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Viajb)["odih"] * (*Rabij)["hfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Viajb)["ofjh"] * (*Rabij)["hdik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Viajb)["ofjh"] * (*Rabij)["heik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Viajb)["oejh"] * (*Rabij)["hdik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Viajb)["odjh"] * (*Rabij)["heik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Viajb)["oejh"] * (*Rabij)["hfik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Viajb)["odjh"] * (*Rabij)["hfik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Viajb)["ofkh"] * (*Rabij)["hdji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Viajb)["ofkh"] * (*Rabij)["heji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Viajb)["oekh"] * (*Rabij)["hdji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Viajb)["odkh"] * (*Rabij)["heji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Viajb)["oekh"] * (*Rabij)["hfji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Viajb)["odkh"] * (*Rabij)["hfji"];

  ST_DEBUG("T1 * Vpppp * R2")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Vabcd)["efhg"] * (*Rabij)["hdjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Vabcd)["fdhg"] * (*Rabij)["hejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Vabcd)["dehg"] * (*Rabij)["hfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Vabcd)["efhg"] * (*Rabij)["hdik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Vabcd)["fdhg"] * (*Rabij)["heik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Vabcd)["dehg"] * (*Rabij)["hfik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Vabcd)["efhg"] * (*Rabij)["hdji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Vabcd)["fdhg"] * (*Rabij)["heji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Vabcd)["dehg"] * (*Rabij)["hfji"];

  ST_DEBUG("T1 * Vhhhp * R3")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Vijka)["pIig"] * (*Rabcijk)["defIpk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Vijka)["pIig"] * (*Rabcijk)["defIpj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Vijka)["pIjg"] * (*Rabcijk)["defIpk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Vijka)["pIkg"] * (*Rabcijk)["defIpj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Vijka)["pIjg"] * (*Rabcijk)["defIpi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Vijka)["pIkg"] * (*Rabcijk)["defIpi"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Vijka)["poiA"] * (*Rabcijk)["Adepjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Vijka)["poiA"] * (*Rabcijk)["Adfpjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Vijka)["poiA"] * (*Rabcijk)["Afepjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Vijka)["pojA"] * (*Rabcijk)["Adepik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Vijka)["pojA"] * (*Rabcijk)["Adfpik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Vijka)["pojA"] * (*Rabcijk)["Afepik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Vijka)["pokA"] * (*Rabcijk)["Adepji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Vijka)["pokA"] * (*Rabcijk)["Adfpji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Vijka)["pokA"] * (*Rabcijk)["Afepji"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gp"] * (*Vijka)["Ipig"] * (*Rabcijk)["defIjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gp"] * (*Vijka)["Ipjg"] * (*Rabcijk)["defIik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gp"] * (*Vijka)["Ipkg"] * (*Rabcijk)["defIji"];

  ST_DEBUG("T1 * Vhppp * R3")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gi"] * (*Viabc)["pfAg"] * (*Rabcijk)["Adepjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Viabc)["peAg"] * (*Rabcijk)["Adfpjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gi"] * (*Viabc)["pdAg"] * (*Rabcijk)["Afepjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gj"] * (*Viabc)["pfAg"] * (*Rabcijk)["Adepik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Viabc)["peAg"] * (*Rabcijk)["Adfpik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gj"] * (*Viabc)["pdAg"] * (*Rabcijk)["Afepik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gk"] * (*Viabc)["pfAg"] * (*Rabcijk)["Adepji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Viabc)["peAg"] * (*Rabcijk)["Adfpji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gk"] * (*Viabc)["pdAg"] * (*Rabcijk)["Afepji"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Viabc)["ofhA"] * (*Rabcijk)["Ahdijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Viabc)["ofhA"] * (*Rabcijk)["Aheijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Viabc)["oehA"] * (*Rabcijk)["Ahdijk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Viabc)["odhA"] * (*Rabcijk)["Aheijk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Viabc)["oehA"] * (*Rabcijk)["Ahfijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Viabc)["odhA"] * (*Rabcijk)["Ahfijk"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gp"] * (*Viabc)["pfAg"] * (*Rabcijk)["Adeijk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gp"] * (*Viabc)["peAg"] * (*Rabcijk)["Adfijk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gp"] * (*Viabc)["pdAg"] * (*Rabcijk)["Afeijk"];

  ST_DEBUG("T2 * Vhhhh * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Vijkl)["poij"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Vijkl)["poij"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Vijkl)["poij"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Vijkl)["poik"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Vijkl)["poik"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Vijkl)["poik"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Vijkl)["pokj"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Vijkl)["pokj"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Vijkl)["pokj"] * (*Rai)["fp"];

  ST_DEBUG("T2 * Vhphp * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gejk"] * (*Viajb)["pfig"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdjk"] * (*Viajb)["pfig"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfjk"] * (*Viajb)["peig"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfjk"] * (*Viajb)["pdig"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdjk"] * (*Viajb)["peig"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gejk"] * (*Viajb)["pdig"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geki"] * (*Viajb)["pfjg"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdki"] * (*Viajb)["pfjg"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfki"] * (*Viajb)["pejg"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfki"] * (*Viajb)["pdjg"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdki"] * (*Viajb)["pejg"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geki"] * (*Viajb)["pdjg"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Viajb)["pfkg"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Viajb)["pfkg"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Viajb)["pekg"] * (*Rai)["dp"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Viajb)["pdkg"] * (*Rai)["ep"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Viajb)["pekg"] * (*Rai)["fp"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Viajb)["pdkg"] * (*Rai)["fp"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Viajb)["ofih"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Viajb)["oeih"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Viajb)["odih"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Viajb)["ofih"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Viajb)["oeih"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Viajb)["odih"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Viajb)["ofjh"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoi"] * (*Viajb)["oejh"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoi"] * (*Viajb)["odjh"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Viajb)["ofkh"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Viajb)["oekh"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Viajb)["odkh"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Viajb)["ofjh"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfok"] * (*Viajb)["oejh"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feok"] * (*Viajb)["odjh"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Viajb)["ofkh"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoj"] * (*Viajb)["oekh"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoj"] * (*Viajb)["odkh"] * (*Rai)["hi"];

  ST_DEBUG("T2 * Vpppp * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Vabcd)["efhg"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Vabcd)["dfhg"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Vabcd)["edhg"] * (*Rai)["hk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Vabcd)["efhg"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Vabcd)["dfhg"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Vabcd)["edhg"] * (*Rai)["hj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Vabcd)["efhg"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Vabcd)["dfhg"] * (*Rai)["hi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Vabcd)["edhg"] * (*Rai)["hi"];

  ST_DEBUG("T2 * Vhhhp * R2")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfjk"] * (*Vijka)["pIig"] * (*Rabij)["deIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gejk"] * (*Vijka)["pIig"] * (*Rabij)["dfIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdjk"] * (*Vijka)["pIig"] * (*Rabij)["feIp"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfki"] * (*Vijka)["pIjg"] * (*Rabij)["deIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geki"] * (*Vijka)["pIjg"] * (*Rabij)["dfIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdki"] * (*Vijka)["pIjg"] * (*Rabij)["feIp"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfij"] * (*Vijka)["pIkg"] * (*Rabij)["deIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geij"] * (*Vijka)["pIkg"] * (*Rabij)["dfIp"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdij"] * (*Vijka)["pIkg"] * (*Rabij)["feIp"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Vijka)["poiA"] * (*Rabij)["Adpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Vijka)["poiA"] * (*Rabij)["Aepk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Vijka)["poiA"] * (*Rabij)["Afpk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Vijka)["poiA"] * (*Rabij)["Adpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Vijka)["poiA"] * (*Rabij)["Aepj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Vijka)["poiA"] * (*Rabij)["Afpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Vijka)["pojA"] * (*Rabij)["Adpk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Vijka)["pojA"] * (*Rabij)["Aepk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Vijka)["pojA"] * (*Rabij)["Afpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Vijka)["pokA"] * (*Rabij)["Adpj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Vijka)["pokA"] * (*Rabij)["Aepj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Vijka)["pokA"] * (*Rabij)["Afpj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Vijka)["pojA"] * (*Rabij)["Adpi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Vijka)["pojA"] * (*Rabij)["Aepi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Vijka)["pojA"] * (*Rabij)["Afpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Vijka)["pokA"] * (*Rabij)["Adpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Vijka)["pokA"] * (*Rabij)["Aepi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Vijka)["pokA"] * (*Rabij)["Afpi"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpj"] * (*Vijka)["Ipig"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepj"] * (*Vijka)["Ipig"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Vijka)["Ipig"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpk"] * (*Vijka)["Ipig"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepk"] * (*Vijka)["Ipig"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpk"] * (*Vijka)["Ipig"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpi"] * (*Vijka)["Ipjg"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepi"] * (*Vijka)["Ipjg"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Vijka)["Ipjg"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpi"] * (*Vijka)["Ipkg"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepi"] * (*Vijka)["Ipkg"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpi"] * (*Vijka)["Ipkg"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpk"] * (*Vijka)["Ipjg"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepk"] * (*Vijka)["Ipjg"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpk"] * (*Vijka)["Ipjg"] * (*Rabij)["feIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpj"] * (*Vijka)["Ipkg"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepj"] * (*Vijka)["Ipkg"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpj"] * (*Vijka)["Ipkg"] * (*Rabij)["feIi"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijka)["opiA"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdop"] * (*Vijka)["opiA"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deop"] * (*Vijka)["opiA"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijka)["opjA"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdop"] * (*Vijka)["opjA"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deop"] * (*Vijka)["opjA"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijka)["opkA"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdop"] * (*Vijka)["opkA"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deop"] * (*Vijka)["opkA"] * (*Rabij)["Afji"];

  ST_DEBUG("T2 * Vhppp * R2")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Viabc)["pfAg"] * (*Rabij)["Adpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Viabc)["pfAg"] * (*Rabij)["Aepk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Viabc)["peAg"] * (*Rabij)["Adpk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Viabc)["pdAg"] * (*Rabij)["Aepk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Viabc)["peAg"] * (*Rabij)["Afpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Viabc)["pdAg"] * (*Rabij)["Afpk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Viabc)["pfAg"] * (*Rabij)["Adpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Viabc)["pfAg"] * (*Rabij)["Aepj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Viabc)["peAg"] * (*Rabij)["Adpj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Viabc)["pdAg"] * (*Rabij)["Aepj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Viabc)["peAg"] * (*Rabij)["Afpj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Viabc)["pdAg"] * (*Rabij)["Afpj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Viabc)["pfAg"] * (*Rabij)["Adpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Viabc)["pfAg"] * (*Rabij)["Aepi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Viabc)["peAg"] * (*Rabij)["Adpi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Viabc)["pdAg"] * (*Rabij)["Aepi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Viabc)["peAg"] * (*Rabij)["Afpi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Viabc)["pdAg"] * (*Rabij)["Afpi"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghij"] * (*Viabc)["Ifgh"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghij"] * (*Viabc)["Iegh"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghij"] * (*Viabc)["Idgh"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghik"] * (*Viabc)["Ifgh"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghik"] * (*Viabc)["Iegh"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghik"] * (*Viabc)["Idgh"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghkj"] * (*Viabc)["Ifgh"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghkj"] * (*Viabc)["Iegh"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghkj"] * (*Viabc)["Idgh"] * (*Rabij)["feIi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deoi"] * (*Viabc)["ofhA"] * (*Rabij)["Ahjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["dfoi"] * (*Viabc)["oehA"] * (*Rabij)["Ahjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["feoi"] * (*Viabc)["odhA"] * (*Rabij)["Ahjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deoj"] * (*Viabc)["ofhA"] * (*Rabij)["Ahik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["dfoj"] * (*Viabc)["oehA"] * (*Rabij)["Ahik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["feoj"] * (*Viabc)["odhA"] * (*Rabij)["Ahik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deok"] * (*Viabc)["ofhA"] * (*Rabij)["Ahji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["dfok"] * (*Viabc)["oehA"] * (*Rabij)["Ahji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["feok"] * (*Viabc)["odhA"] * (*Rabij)["Ahji"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepi"] * (*Viabc)["pfAg"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Viabc)["pfAg"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpi"] * (*Viabc)["peAg"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpi"] * (*Viabc)["pdAg"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpi"] * (*Viabc)["peAg"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepi"] * (*Viabc)["pdAg"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepj"] * (*Viabc)["pfAg"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Viabc)["pfAg"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpj"] * (*Viabc)["peAg"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpj"] * (*Viabc)["pdAg"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpj"] * (*Viabc)["peAg"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepj"] * (*Viabc)["pdAg"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepk"] * (*Viabc)["pfAg"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpk"] * (*Viabc)["pfAg"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpk"] * (*Viabc)["peAg"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpk"] * (*Viabc)["pdAg"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpk"] * (*Viabc)["peAg"] * (*Rabij)["Afji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepk"] * (*Viabc)["pdAg"] * (*Rabij)["Afji"];

  ST_DEBUG("T2 * Vhhpp * R3")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfij"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdeIpk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geij"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdfIpk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdij"] * (*Vijab)["pIBg"] * (*Rabcijk)["BfeIpk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gfik"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdeIpj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["geik"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdfIpj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gdik"] * (*Vijab)["pIBg"] * (*Rabcijk)["BfeIpj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gfkj"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdeIpi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gekj"] * (*Vijab)["pIBg"] * (*Rabcijk)["BdfIpi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gdkj"] * (*Vijab)["pIBg"] * (*Rabcijk)["BfeIpi"];

  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabij)["ghij"] * (*Vijab)["IJgh"] * (*Rabcijk)["defJIk"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabij)["ghik"] * (*Vijab)["IJgh"] * (*Rabcijk)["defJIj"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabij)["ghkj"] * (*Vijab)["IJgh"] * (*Rabcijk)["defJIi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efoi"] * (*Vijab)["poAB"] * (*Rabcijk)["BAdpjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdoi"] * (*Vijab)["poAB"] * (*Rabcijk)["BAepjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deoi"] * (*Vijab)["poAB"] * (*Rabcijk)["BAfpjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efoj"] * (*Vijab)["poAB"] * (*Rabcijk)["BAdpik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdoj"] * (*Vijab)["poAB"] * (*Rabcijk)["BAepik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deoj"] * (*Vijab)["poAB"] * (*Rabcijk)["BAfpik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efok"] * (*Vijab)["poAB"] * (*Rabcijk)["BAdpji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdok"] * (*Vijab)["poAB"] * (*Rabcijk)["BAepji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deok"] * (*Vijab)["poAB"] * (*Rabcijk)["BAfpji"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpi"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdeIjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepi"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdfIjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpi"] * (*Vijab)["IpBg"] * (*Rabcijk)["BfeIjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpj"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdeIik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepj"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdfIik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpj"] * (*Vijab)["IpBg"] * (*Rabcijk)["BfeIik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpk"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdeIji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepk"] * (*Vijab)["IpBg"] * (*Rabcijk)["BdfIji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpk"] * (*Vijab)["IpBg"] * (*Rabcijk)["BfeIji"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghIi"] * (*Vijab)["JIgh"] * (*Rabcijk)["defJjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghIj"] * (*Vijab)["JIgh"] * (*Rabcijk)["defJik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghIk"] * (*Vijab)["JIgh"] * (*Rabcijk)["defJji"];

  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabij)["efop"] * (*Vijab)["opAB"] * (*Rabcijk)["BAdijk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabij)["fdop"] * (*Vijab)["opAB"] * (*Rabcijk)["BAeijk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabij)["deop"] * (*Vijab)["opAB"] * (*Rabcijk)["BAfijk"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gfpI"] * (*Vijab)["pIBg"] * (*Rabcijk)["Bdeijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gepI"] * (*Vijab)["pIBg"] * (*Rabcijk)["Bdfijk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gdpI"] * (*Vijab)["pIBg"] * (*Rabcijk)["Bfeijk"];

  ST_DEBUG("T1 * T1 * Vhhhp * R2")
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rabij)["feIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rabij)["feIi"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Vijka)["poiA"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["poiA"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Vijka)["poiA"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["poiA"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Vijka)["poiA"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["poiA"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Vijka)["pojA"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["pojA"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Vijka)["pojA"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["pojA"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Vijka)["pojA"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["pojA"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Vijka)["pokA"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijka)["pokA"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Vijka)["pokA"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Vijka)["pokA"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Vijka)["pokA"] * (*Rabij)["Afji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Vijka)["pokA"] * (*Rabij)["Afji"];

  ST_DEBUG("T1 * T1 * Vhppp * R2")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Ifhg"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Viabc)["Ifhg"] * (*Rabij)["deIk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Iehg"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Viabc)["Iehg"] * (*Rabij)["dfIk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Viabc)["Idhg"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Viabc)["Idhg"] * (*Rabij)["feIk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Ifhg"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Viabc)["Ifhg"] * (*Rabij)["deIj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Iehg"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Viabc)["Iehg"] * (*Rabij)["dfIj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Viabc)["Idhg"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Viabc)["Idhg"] * (*Rabij)["feIj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Ifhg"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Viabc)["Ifhg"] * (*Rabij)["deIi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Iehg"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Viabc)["Iehg"] * (*Rabij)["dfIi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Viabc)["Idhg"] * (*Rabij)["feIi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Viabc)["Idhg"] * (*Rabij)["feIi"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Viabc)["ofAh"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Viabc)["ofAh"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Viabc)["oeAh"] * (*Rabij)["Adjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Viabc)["odAh"] * (*Rabij)["Aejk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Viabc)["oeAh"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Viabc)["odAh"] * (*Rabij)["Afjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Viabc)["ofAh"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Viabc)["ofAh"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Viabc)["oeAh"] * (*Rabij)["Adik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Viabc)["odAh"] * (*Rabij)["Aeik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Viabc)["oeAh"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Viabc)["odAh"] * (*Rabij)["Afik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Viabc)["ofAh"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Viabc)["ofAh"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Viabc)["oeAh"] * (*Rabij)["Adji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Viabc)["odAh"] * (*Rabij)["Aeji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Viabc)["oeAh"] * (*Rabij)["Afji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Viabc)["odAh"] * (*Rabij)["Afji"];

  ST_DEBUG("T1 * T1 * Vhhpp * R3")
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIj"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIj"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIi"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Vijab)["IJhg"] * (*Rabcijk)["defJIi"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdeIjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdfIjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabcijk)["BfeIjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdeIik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdfIik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabcijk)["BfeIik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdeIji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabcijk)["BdfIji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabcijk)["BfeIji"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["gi"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabcijk)["defJjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gj"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabcijk)["defJik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["gk"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabcijk)["defJji"];

  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Vijab)["poAB"] * (*Rabcijk)["BAdijk"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Vijab)["poAB"] * (*Rabcijk)["BAdijk"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Vijab)["poAB"] * (*Rabcijk)["BAeijk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Vijab)["poAB"] * (*Rabcijk)["BAeijk"];
  (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Vijab)["poAB"] * (*Rabcijk)["BAfijk"];
  (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Vijab)["poAB"] * (*Rabcijk)["BAfijk"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabcijk)["Bdeijk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabcijk)["Bdfijk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabcijk)["Bfeijk"];

  ST_DEBUG("T2 * T1 * Vhhhp * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Vijka)["Ioih"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Vijka)["Ioih"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Vijka)["Iojh"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Vijka)["Iokh"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Vijka)["Iojh"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Vijka)["Iokh"] * (*Rai)["fI"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gejk"] * (*Tai)["fp"] * (*Vijka)["Ipig"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["fp"] * (*Vijka)["Ipig"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["ep"] * (*Vijka)["Ipig"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["dp"] * (*Vijka)["Ipig"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["ep"] * (*Vijka)["Ipig"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["dp"] * (*Vijka)["Ipig"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geki"] * (*Tai)["fp"] * (*Vijka)["Ipjg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdki"] * (*Tai)["fp"] * (*Vijka)["Ipjg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["ep"] * (*Vijka)["Ipjg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfki"] * (*Tai)["dp"] * (*Vijka)["Ipjg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["ep"] * (*Vijka)["Ipjg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["dp"] * (*Vijka)["Ipjg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Vijka)["Ipkg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Vijka)["Ipkg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Vijka)["Ipkg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Vijka)["Ipkg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Vijka)["Ipkg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Vijka)["Ipkg"] * (*Rai)["fI"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijka)["poiA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijka)["poiA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijka)["poiA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijka)["poiA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijka)["poiA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijka)["poiA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijka)["pojA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijka)["pojA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijka)["pojA"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijka)["pokA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijka)["pokA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijka)["pokA"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijka)["pojA"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijka)["pojA"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijka)["pojA"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijka)["pokA"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijka)["pokA"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijka)["pokA"] * (*Rai)["Ai"];

  ST_DEBUG("T2 * T1 * Vhppp * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Viabc)["Ifhg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Viabc)["Ifhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Viabc)["Iehg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Viabc)["Idhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Viabc)["Iehg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Viabc)["Idhg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Viabc)["Ifhg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Viabc)["Ifhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Viabc)["Iehg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Viabc)["Idhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Viabc)["Iehg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Viabc)["Idhg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Viabc)["Ifhg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Viabc)["Ifhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Viabc)["Iehg"] * (*Rai)["dI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Viabc)["Idhg"] * (*Rai)["eI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Viabc)["Iehg"] * (*Rai)["fI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Viabc)["Idhg"] * (*Rai)["fI"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Viabc)["ofAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hi"] * (*Viabc)["oeAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hi"] * (*Viabc)["odAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Viabc)["ofAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hi"] * (*Viabc)["oeAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feok"] * (*Tai)["hi"] * (*Viabc)["odAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Viabc)["ofAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hj"] * (*Viabc)["oeAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hj"] * (*Viabc)["odAh"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Viabc)["ofAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfoi"] * (*Tai)["hk"] * (*Viabc)["oeAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feoi"] * (*Tai)["hk"] * (*Viabc)["odAh"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Viabc)["ofAh"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Tai)["hj"] * (*Viabc)["oeAh"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Tai)["hj"] * (*Viabc)["odAh"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Viabc)["ofAh"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Tai)["hk"] * (*Viabc)["oeAh"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Tai)["hk"] * (*Viabc)["odAh"] * (*Rai)["Ai"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Viabc)["pfAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Viabc)["pfAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Viabc)["peAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Viabc)["pdAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Viabc)["peAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Viabc)["pdAg"] * (*Rai)["Ak"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tai)["ep"] * (*Viabc)["pfAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tai)["dp"] * (*Viabc)["pfAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tai)["fp"] * (*Viabc)["peAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tai)["fp"] * (*Viabc)["pdAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tai)["dp"] * (*Viabc)["peAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tai)["ep"] * (*Viabc)["pdAg"] * (*Rai)["Aj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tai)["ep"] * (*Viabc)["pfAg"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tai)["dp"] * (*Viabc)["pfAg"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tai)["fp"] * (*Viabc)["peAg"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tai)["fp"] * (*Viabc)["pdAg"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tai)["dp"] * (*Viabc)["peAg"] * (*Rai)["Ai"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tai)["ep"] * (*Viabc)["pdAg"] * (*Rai)["Ai"];

  ST_DEBUG("T2 * T1 * Vhhpp * R2")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Vijab)["IJhg"] * (*Rabij)["deJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Vijab)["IJhg"] * (*Rabij)["dfJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Vijab)["IJhg"] * (*Rabij)["feJI"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Vijab)["IJhg"] * (*Rabij)["deJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Vijab)["IJhg"] * (*Rabij)["dfJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Vijab)["IJhg"] * (*Rabij)["feJI"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Vijab)["IJhg"] * (*Rabij)["deJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Vijab)["IJhg"] * (*Rabij)["dfJI"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Vijab)["IJhg"] * (*Rabij)["feJI"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BdIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BeIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BfIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BdIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BeIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rabij)["BfIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BdIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BeIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BfIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BdIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BeIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BfIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BdIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BeIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rabij)["BfIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BdIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BeIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rabij)["BfIi"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpj"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepj"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpi"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepi"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpi"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepi"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpk"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepk"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpk"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rabij)["feJi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpj"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepj"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rabij)["feJi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Tai)["Ai"] * (*Vijab)["opBA"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdop"] * (*Tai)["Ai"] * (*Vijab)["opBA"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deop"] * (*Tai)["Ai"] * (*Vijab)["opBA"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Tai)["Aj"] * (*Vijab)["opBA"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdop"] * (*Tai)["Aj"] * (*Vijab)["opBA"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deop"] * (*Tai)["Aj"] * (*Vijab)["opBA"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Tai)["Ak"] * (*Vijab)["opBA"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdop"] * (*Tai)["Ak"] * (*Vijab)["opBA"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deop"] * (*Tai)["Ak"] * (*Vijab)["opBA"] * (*Rabij)["Bfji"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BdIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BdIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BfIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BfIk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BdIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BdIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BfIj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BfIj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BdIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BdIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BeIi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rabij)["BfIi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rabij)["BfIi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghij"] * (*Tai)["fI"] * (*Vijab)["JIgh"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghij"] * (*Tai)["eI"] * (*Vijab)["JIgh"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghij"] * (*Tai)["dI"] * (*Vijab)["JIgh"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghik"] * (*Tai)["fI"] * (*Vijab)["JIgh"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghik"] * (*Tai)["eI"] * (*Vijab)["JIgh"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghik"] * (*Tai)["dI"] * (*Vijab)["JIgh"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["ghkj"] * (*Tai)["fI"] * (*Vijab)["JIgh"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghkj"] * (*Tai)["eI"] * (*Vijab)["JIgh"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["ghkj"] * (*Tai)["dI"] * (*Vijab)["JIgh"] * (*Rabij)["feJi"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rabij)["feJi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijab)["poAB"] * (*Rabij)["BAjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijab)["poAB"] * (*Rabij)["BAjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijab)["poAB"] * (*Rabij)["BAjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijab)["poAB"] * (*Rabij)["BAik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijab)["poAB"] * (*Rabij)["BAik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijab)["poAB"] * (*Rabij)["BAik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijab)["poAB"] * (*Rabij)["BAji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijab)["poAB"] * (*Rabij)["BAji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijab)["poAB"] * (*Rabij)["BAji"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepi"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpi"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpi"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepi"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepj"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpj"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpj"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepj"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gepk"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdpk"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfpk"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfpk"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdpk"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfji"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gepk"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rabij)["Bfji"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rabij)["Bfji"];

  ST_DEBUG("T2 * T2 * Vhhpp * R1")
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfij"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["dfok"] * (*Tabij)["heij"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["feok"] * (*Tabij)["hdij"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfik"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoj"] * (*Tabij)["heik"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoj"] * (*Tabij)["hdik"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["feoi"] * (*Tabij)["hdkj"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["dfoi"] * (*Tabij)["hekj"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfkj"] * (*Vijab)["IoBh"] * (*Rai)["BI"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["efok"] * (*Tabij)["hAij"] * (*Vijab)["JohA"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["fdok"] * (*Tabij)["hAij"] * (*Vijab)["JohA"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Vijab)["JohA"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efoj"] * (*Tabij)["hAik"] * (*Vijab)["JohA"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdoj"] * (*Tabij)["hAik"] * (*Vijab)["JohA"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Vijab)["JohA"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["efoi"] * (*Tabij)["hAkj"] * (*Vijab)["JohA"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["fdoi"] * (*Tabij)["hAkj"] * (*Vijab)["JohA"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Vijab)["JohA"] * (*Rai)["fJ"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Vijab)["pIBg"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Vijab)["pIBg"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Vijab)["pIBg"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfik"] * (*Tabij)["depI"] * (*Vijab)["pIBg"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["geik"] * (*Tabij)["dfpI"] * (*Vijab)["pIBg"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdik"] * (*Tabij)["fepI"] * (*Vijab)["pIBg"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabij)["gfkj"] * (*Tabij)["depI"] * (*Vijab)["pIBg"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gekj"] * (*Tabij)["dfpI"] * (*Vijab)["pIBg"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabij)["gdkj"] * (*Tabij)["fepI"] * (*Vijab)["pIBg"] * (*Rai)["Bi"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];

  ST_DEBUG("T1 * T1 * T1 * Vhhpp * R2")
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["deJk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["dfJk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["feJk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["deJj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["dfJj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tai)["Ai"] * (*Vijab)["JoAh"] * (*Rabij)["feJj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["deJi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["dfJi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tai)["Ak"] * (*Vijab)["JoAh"] * (*Rabij)["feJi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tai)["Aj"] * (*Vijab)["JoAh"] * (*Rabij)["feJi"];

  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bdjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bejk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tai)["Ai"] * (*Vijab)["poBA"] * (*Rabij)["Bfjk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Bdik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Beik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tai)["Aj"] * (*Vijab)["poBA"] * (*Rabij)["Bfik"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Bdji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Beji"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Bfji"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tai)["Ak"] * (*Vijab)["poBA"] * (*Rabij)["Bfji"];

  ST_DEBUG("T1 * T1 * T2 * V * R1")
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["efIk"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Tabij)["efIk"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["fdIk"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Tabij)["fdIk"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tabij)["deIk"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gj"] * (*Tai)["hi"] * (*Tabij)["deIk"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["efIj"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Tabij)["efIj"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["fdIj"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Tabij)["fdIj"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tabij)["deIj"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gk"] * (*Tai)["hi"] * (*Tabij)["deIj"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["efIi"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Tabij)["efIi"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["fdIi"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Tabij)["fdIi"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tabij)["deIi"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["gk"] * (*Tai)["hj"] * (*Tabij)["deIi"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];

  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tabij)["Aejk"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tabij)["Adjk"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tabij)["Afjk"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tabij)["Afjk"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tabij)["Adjk"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tabij)["Aejk"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tabij)["Aeik"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tabij)["Adik"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tabij)["Afik"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tabij)["Afik"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tabij)["Adik"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tabij)["Aeik"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tabij)["Aeji"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tabij)["Adji"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tabij)["Afji"] * (*Vijab)["JoAh"] * (*Rai)["dJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tabij)["Afji"] * (*Vijab)["JoAh"] * (*Rai)["eJ"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tabij)["Adji"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tabij)["Aeji"] * (*Vijab)["JoAh"] * (*Rai)["fJ"];

  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tabij)["deIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tabij)["dfIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tabij)["feIj"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hi"] * (*Tabij)["deIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hi"] * (*Tabij)["dfIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hi"] * (*Tabij)["feIk"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tabij)["deIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tabij)["dfIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tabij)["feIi"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tabij)["deIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tabij)["dfIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tabij)["feIi"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["fo"] * (*Tai)["hj"] * (*Tabij)["deIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["eo"] * (*Tai)["hj"] * (*Tabij)["dfIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["do"] * (*Tai)["hj"] * (*Tabij)["feIk"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tai)["fo"] * (*Tai)["hk"] * (*Tabij)["deIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["eo"] * (*Tai)["hk"] * (*Tabij)["dfIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tai)["do"] * (*Tai)["hk"] * (*Tabij)["feIj"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];

  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tabij)["Adij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tabij)["Aeij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aeij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tabij)["Afij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afij"] * (*Vijab)["poBA"] * (*Rai)["Bk"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tabij)["Adik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tabij)["Aeik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aeik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tabij)["Afik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afik"] * (*Vijab)["poBA"] * (*Rai)["Bj"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["fo"] * (*Tai)["ep"] * (*Tabij)["Adkj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["eo"] * (*Tai)["fp"] * (*Tabij)["Adkj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["fo"] * (*Tai)["dp"] * (*Tabij)["Aekj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["do"] * (*Tai)["fp"] * (*Tabij)["Aekj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tai)["eo"] * (*Tai)["dp"] * (*Tabij)["Afkj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];
  (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tai)["do"] * (*Tai)["ep"] * (*Tabij)["Afkj"] * (*Vijab)["poBA"] * (*Rai)["Bi"];

  if (dressing == CCSDT) {
    ST_DEBUG("T3 * * V * R1")
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["defojk"] * (*Vijka)["poiA"] * (*Rai)["Ap"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["defoki"] * (*Vijka)["pojA"] * (*Rai)["Ap"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["defoij"] * (*Vijka)["pokA"] * (*Rai)["Ap"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpjk"] * (*Vijka)["Ipig"] * (*Rai)["dI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpjk"] * (*Vijka)["Ipig"] * (*Rai)["eI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepjk"] * (*Vijka)["Ipig"] * (*Rai)["fI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpki"] * (*Vijka)["Ipjg"] * (*Rai)["dI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpki"] * (*Vijka)["Ipjg"] * (*Rai)["eI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepki"] * (*Vijka)["Ipjg"] * (*Rai)["fI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpij"] * (*Vijka)["Ipkg"] * (*Rai)["dI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpij"] * (*Vijka)["Ipkg"] * (*Rai)["eI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Vijka)["Ipkg"] * (*Rai)["fI"];

    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Vijka)["opiA"] * (*Rai)["Ak"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Vijka)["opiA"] * (*Rai)["Aj"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopi"] * (*Vijka)["opjA"] * (*Rai)["Ak"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Vijka)["opkA"] * (*Rai)["Aj"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopk"] * (*Vijka)["opjA"] * (*Rai)["Ai"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopj"] * (*Vijka)["opkA"] * (*Rai)["Ai"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdeijk"] * (*Viabc)["pfAg"] * (*Rai)["Ap"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdfijk"] * (*Viabc)["peAg"] * (*Rai)["Ap"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfeijk"] * (*Viabc)["pdAg"] * (*Rai)["Ap"];

    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gheijk"] * (*Viabc)["Ifgh"] * (*Rai)["dI"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghdijk"] * (*Viabc)["Ifgh"] * (*Rai)["eI"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Viabc)["Iegh"] * (*Rai)["dI"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghfijk"] * (*Viabc)["Idgh"] * (*Rai)["eI"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Viabc)["Iegh"] * (*Rai)["fI"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Viabc)["Idgh"] * (*Rai)["fI"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Viabc)["pfAg"] * (*Rai)["Ak"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdfpij"] * (*Viabc)["peAg"] * (*Rai)["Ak"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gfepij"] * (*Viabc)["pdAg"] * (*Rai)["Ak"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdepik"] * (*Viabc)["pfAg"] * (*Rai)["Aj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpik"] * (*Viabc)["peAg"] * (*Rai)["Aj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfepik"] * (*Viabc)["pdAg"] * (*Rai)["Aj"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdepkj"] * (*Viabc)["pfAg"] * (*Rai)["Ai"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpkj"] * (*Viabc)["peAg"] * (*Rai)["Ai"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfepkj"] * (*Viabc)["pdAg"] * (*Rai)["Ai"];

    ST_DEBUG("T3 * V * R2")
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gefijk"] * (*Vijab)["pIBg"] * (*Rabij)["BdIp"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gfdijk"] * (*Vijab)["pIBg"] * (*Rabij)["BeIp"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gdeijk"] * (*Vijab)["pIBg"] * (*Rabij)["BfIp"];

    (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabcijk)["ghfijk"] * (*Vijab)["IJgh"] * (*Rabij)["deJI"];
    (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabcijk)["gheijk"] * (*Vijab)["IJgh"] * (*Rabij)["dfJI"];
    (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabcijk)["ghdijk"] * (*Vijab)["IJgh"] * (*Rabij)["feJI"];

    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defoij"] * (*Vijab)["poAB"] * (*Rabij)["BApk"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defoik"] * (*Vijab)["poAB"] * (*Rabij)["BApj"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defokj"] * (*Vijab)["poAB"] * (*Rabij)["BApi"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gefpij"] * (*Vijab)["IpBg"] * (*Rabij)["BdIk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gfdpij"] * (*Vijab)["IpBg"] * (*Rabij)["BeIk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdepij"] * (*Vijab)["IpBg"] * (*Rabij)["BfIk"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpik"] * (*Vijab)["IpBg"] * (*Rabij)["BdIj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpik"] * (*Vijab)["IpBg"] * (*Rabij)["BeIj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepik"] * (*Vijab)["IpBg"] * (*Rabij)["BfIj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpkj"] * (*Vijab)["IpBg"] * (*Rabij)["BdIi"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpkj"] * (*Vijab)["IpBg"] * (*Rabij)["BeIi"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepkj"] * (*Vijab)["IpBg"] * (*Rabij)["BfIi"];

    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghfIij"] * (*Vijab)["JIgh"] * (*Rabij)["deJk"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gheIij"] * (*Vijab)["JIgh"] * (*Rabij)["dfJk"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghdIij"] * (*Vijab)["JIgh"] * (*Rabij)["feJk"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghfIik"] * (*Vijab)["JIgh"] * (*Rabij)["deJj"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gheIik"] * (*Vijab)["JIgh"] * (*Rabij)["dfJj"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghdIik"] * (*Vijab)["JIgh"] * (*Rabij)["feJj"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghfIkj"] * (*Vijab)["JIgh"] * (*Rabij)["deJi"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gheIkj"] * (*Vijab)["JIgh"] * (*Rabij)["dfJi"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghdIkj"] * (*Vijab)["JIgh"] * (*Rabij)["feJi"];

    (*HRabcijk)["defijk"] += ( - 0.25  ) * (*Tabcijk)["defopi"] * (*Vijab)["opAB"] * (*Rabij)["BAjk"];
    (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabcijk)["defopj"] * (*Vijab)["opAB"] * (*Rabij)["BAik"];
    (*HRabcijk)["defijk"] += ( + 0.25  ) * (*Tabcijk)["defopk"] * (*Vijab)["opAB"] * (*Rabij)["BAji"];

    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gefpIi"] * (*Vijab)["pIBg"] * (*Rabij)["Bdjk"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gfdpIi"] * (*Vijab)["pIBg"] * (*Rabij)["Bejk"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gdepIi"] * (*Vijab)["pIBg"] * (*Rabij)["Bfjk"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gefpIj"] * (*Vijab)["pIBg"] * (*Rabij)["Bdik"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gfdpIj"] * (*Vijab)["pIBg"] * (*Rabij)["Beik"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gdepIj"] * (*Vijab)["pIBg"] * (*Rabij)["Bfik"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gefpIk"] * (*Vijab)["pIBg"] * (*Rabij)["Bdji"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gfdpIk"] * (*Vijab)["pIBg"] * (*Rabij)["Beji"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gdepIk"] * (*Vijab)["pIBg"] * (*Rabij)["Bfji"];

    ST_DEBUG("T3 * T1 * V * R1")
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["defojk"] * (*Tai)["hi"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["defoki"] * (*Tai)["hj"] * (*Vijab)["IoBh"] * (*Rai)["BI"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["defoij"] * (*Tai)["hk"] * (*Vijab)["IoBh"] * (*Rai)["BI"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpjk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpjk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepjk"] * (*Tai)["Ai"] * (*Vijab)["JpAg"] * (*Rai)["fJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpki"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpki"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepki"] * (*Tai)["Aj"] * (*Vijab)["JpAg"] * (*Rai)["fJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gefpij"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfdpij"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["Ak"] * (*Vijab)["JpAg"] * (*Rai)["fJ"];

    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ai"] * (*Vijab)["opBA"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Ai"] * (*Vijab)["opBA"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Aj"] * (*Vijab)["opBA"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopi"] * (*Tai)["Ak"] * (*Vijab)["opBA"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["defopk"] * (*Tai)["Aj"] * (*Vijab)["opBA"] * (*Rai)["Bi"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["defopj"] * (*Tai)["Ak"] * (*Vijab)["opBA"] * (*Rai)["Bi"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdeijk"] * (*Tai)["fp"] * (*Vijab)["IpBg"] * (*Rai)["BI"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdfijk"] * (*Tai)["ep"] * (*Vijab)["IpBg"] * (*Rai)["BI"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gfeijk"] * (*Tai)["dp"] * (*Vijab)["IpBg"] * (*Rai)["BI"];

    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["fI"] * (*Vijab)["JIgh"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["fI"] * (*Vijab)["JIgh"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["eI"] * (*Vijab)["JIgh"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghfijk"] * (*Tai)["dI"] * (*Vijab)["JIgh"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( - 0.5  ) * (*Tabcijk)["ghdijk"] * (*Tai)["eI"] * (*Vijab)["JIgh"] * (*Rai)["fJ"];
    (*HRabcijk)["defijk"] += ( + 0.5  ) * (*Tabcijk)["gheijk"] * (*Tai)["dI"] * (*Vijab)["JIgh"] * (*Rai)["fJ"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gefijk"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rai)["dJ"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gfdijk"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rai)["eJ"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdeijk"] * (*Tai)["hI"] * (*Vijab)["JIhg"] * (*Rai)["fJ"];

    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdepij"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdfpij"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gfepij"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdepik"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpik"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfepik"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["gdepkj"] * (*Tai)["fI"] * (*Vijab)["IpBg"] * (*Rai)["Bi"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gdfpkj"] * (*Tai)["eI"] * (*Vijab)["IpBg"] * (*Rai)["Bi"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["gfepkj"] * (*Tai)["dI"] * (*Vijab)["IpBg"] * (*Rai)["Bi"];

    (*HRabcijk)["defijk"] += ( + 1.0  ) * (*Tabcijk)["defoij"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rai)["Bk"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["defoik"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rai)["Bj"];
    (*HRabcijk)["defijk"] += ( - 1.0  ) * (*Tabcijk)["defokj"] * (*Tai)["hI"] * (*Vijab)["IoBh"] * (*Rai)["Bi"];
  }

  ST_DEBUG("triples done")

  return HR;

}

template <typename F>
PTR(StantonIntermediatesUCCSD<F>)
SimilarityTransformedHamiltonian<F>::getStantonIntermediatesUCCSD() {
  if (stantonIntermediatesUccsd) return stantonIntermediatesUccsd;

  stantonIntermediatesUccsd = NEW(StantonIntermediatesUCCSD<F>);
  stantonIntermediatesUccsd->setTai(Tai);
  stantonIntermediatesUccsd->setTabij(Tabij);
  stantonIntermediatesUccsd->setVabcd(Vabcd);
  stantonIntermediatesUccsd->setViajb(Viajb);
  stantonIntermediatesUccsd->setVijab(Vijab);
  stantonIntermediatesUccsd->setVijkl(Vijkl);
  stantonIntermediatesUccsd->setVijka(Vijka);
  stantonIntermediatesUccsd->setViabc(Viabc);
  stantonIntermediatesUccsd->setViajk(Viajk);
  stantonIntermediatesUccsd->setVabic(Vabic);
  stantonIntermediatesUccsd->setVaibc(Vaibc);
  stantonIntermediatesUccsd->setVaibj(Vaibj);
  stantonIntermediatesUccsd->setViabj(Viabj);
  stantonIntermediatesUccsd->setVijak(Vijak);
  stantonIntermediatesUccsd->setVaijb(Vaijb);
  stantonIntermediatesUccsd->setVabci(Vabci);
  stantonIntermediatesUccsd->setVabij(Vabij);
  stantonIntermediatesUccsd->setFij(Fij);
  stantonIntermediatesUccsd->setFab(Fab);
  if (Fia) stantonIntermediatesUccsd->setFia(Fia);
  return stantonIntermediatesUccsd;
}

// instantiate
template
class SimilarityTransformedHamiltonian<cc4s::complex>;
template
class SimilarityTransformedHamiltonian<double>;
