#include <algorithms/StantonIntermediatesUCCSD.hpp>
#include <vector>
#include <map>
#include <ctf.hpp>
#include <string>
#include <util/Exception.hpp>
#include <util/SharedPointer.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

template <typename F>
void StantonIntermediatesUCCSD<F>::checkInputs() {

  // Check if everything needed is given
  std::map<std::string, CTF::Tensor<F>*> inputTensors{
    {"Tai", Tai}, {"Tabij", Tabij},
    {"Fij", Fij}, {"Fab", Fab},
    {"Vabcd", Vabcd}, {"Viajb", Viajb}, {"Vijab", Vijab},
    {"Vijkl", Vijkl}, {"Vijka", Vijka}, {"Viabc", Viabc},
    {"Viajk", Viajk}, {"Vabic", Vabic}, {"Vaibc", Vaibc},
    {"Vaibj", Vaibj}, {"Viabj", Viabj}, {"Vijak", Vijak},
    {"Vaijb", Vaijb}, {"Vabci", Vabci}, {"Vabij", Vabij}
  };

  for (auto entry: inputTensors) {
    if (! entry.second) {
      std::string message("You need: ");
      message.append(entry.first);
      throw new EXCEPTION(message);
    }
  }

}

template <typename F>
void StantonIntermediatesUCCSD<F>::calculateOneBodyIntermediates() {

  if (Fae && Fmi && Fme && TildeTau_abij) { return; }

  checkInputs();

  // Define intermediates
  int No(Fij->lens[0]), Nv(Fab->lens[0]);
  int vv[] = {Nv, Nv}, oo[] = {No, No}, ov[] = {No, Nv};
  int syms[] = {NS, NS};

  Fae = NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "Fae");
  Fmi = NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "Fmi");
  Fme = NEW(CTF::Tensor<F>, 2, ov, syms, *Cc4s::world, "Fme");

  // Equation (9)
  TildeTau_abij = NEW(CTF::Tensor<F>, *Tabij);
  (*TildeTau_abij)["abij"] += ( 0.5 ) * (*Tai)["ai"] * (*Tai)["bj"];
  (*TildeTau_abij)["abij"] += ( - 0.5 ) * (*Tai)["bi"] * (*Tai)["aj"];

  // Equation (3)
  (*Fae)["ae"]  = (*Fab)["ae"];
  (*Fae)["aa"] += (-1.0) * (*Fab)["aa"];
  if (Fia) {
    (*Fae)["aa"] += (-0.5) * (*Fia)["me"] * (*Tai)["am"];
  }
  (*Fae)["ae"] += (*Tai)["fm"] * (*Viabc)["mafe"];
  (*Fae)["ae"] += ( - 0.5 ) * (*TildeTau_abij)["afmn"] * (*Vijab)["mnef"];

  // Equation (4)
  (*Fmi)["mi"]  = (*Fij)["mi"];
  (*Fmi)["ii"] += (-1.0) * (*Fij)["ii"];
  if (Fia) {
    (*Fmi)["mi"] += (+0.5) * (*Fia)["me"] * (*Tai)["ei"];
  }
  (*Fmi)["mi"] += (*Tai)["en"] * (*Vijka)["mnie"];
  (*Fmi)["mi"] += (0.5) * (*TildeTau_abij)["efin"] * (*Vijab)["mnef"];

  // Equation (5) Stanton et al.
  (*Fme)["me"] = (*Tai)["fn"] * (*Vijab)["mnef"];
  if (Fia) {
    (*Fme)["me"] += (*Fia)["me"];
  }

}

template <typename F>
void StantonIntermediatesUCCSD<F>::calculateIntermediates() {

  calculateOneBodyIntermediates();
  if (Tau_abij && Wijkl && Wabcd && Wiabj) { return; }
  checkInputs();

  // Equation (10)
  Tau_abij = NEW(CTF::Tensor<F>, *Tabij);
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

  // Equation (6)
  Wijkl = NEW(CTF::Tensor<F>, *Vijkl);
  (*Wijkl)["mnij"] += (+ 1.0) * (*Tai)["ej"] * (*Vijka)["mnie"];
  // Pij
  (*Wijkl)["mnij"] += (- 1.0) * (*Tai)["ei"] * (*Vijka)["mnje"];
  (*Wijkl)["mnij"] += (0.25) * (*Tau_abij)["efij"] * (*Vijab)["mnef"];

  // Equation (7)
  Wabcd = NEW(CTF::Tensor<F>, *Vabcd);
  (*Wabcd)["abef"] += (- 1.0) * (*Tai)["bm"] * (*Vaibc)["amef"];
  // Pab
  (*Wabcd)["abef"] += (+ 1.0) * (*Tai)["am"] * (*Vaibc)["bmef"];
  (*Wabcd)["abef"] += (0.25) * (*Tau_abij)["abmn"] * (*Vijab)["mnef"];

  // Equation (8)
  Wiabj = NEW(CTF::Tensor<F>, *Viabj);
  (*Wiabj)["mbej"] += (+ 1.0) * (*Tai)["fj"] * (*Viabc)["mbef"];
  (*Wiabj)["mbej"] += (- 1.0) * (*Tai)["bn"] * (*Vijak)["mnej"];
  (*Wiabj)["mbej"] += (- 0.5) * (*Tabij)["fbjn"] * (*Vijab)["mnef"];
  (*Wiabj)["mbej"] += (- 1.0) * (*Tai)["fj"] * (*Tai)["bn"] * (*Vijab)["mnef"];

}

template <typename F>
PTR(CTF::Tensor<F>)
StantonIntermediatesUCCSD<F>::getRai(){

  if (Rai) return Rai;

  calculateOneBodyIntermediates();

  int No(Fij->lens[0]), Nv(Fab->lens[0]);
  int vo[] = {Nv, No};
  int syms[] = {NS, NS};

  Rai = NEW(CTF::Tensor<F>, 2, vo, syms, *Cc4s::world, "Rai");

  // T1 equations:
  (*Rai)["ai"] = (*Tai)["ei"] * (*Fae)["ae"];
  if (Fia) {
     (*Rai)["ai"] += (*Fia)["ia"] ;
  }

  (*Rai)["ai"] += (- 1.0) * (*Tai)["am"] * (*Fmi)["mi"];
  (*Rai)["ai"] += (*Tabij)["aeim"] * (*Fme)["me"];
  (*Rai)["ai"] += (- 1.0) * (*Tai)["fn"] * (*Viajb)["naif"];
  (*Rai)["ai"] += (- 0.5) * (*Tabij)["efim"] * (*Viabc)["maef"];
  (*Rai)["ai"] += (- 0.5) * (*Tabij)["aemn"] * (*Vijak)["nmei"];

  return Rai;

}

template <typename F>
PTR(CTF::Tensor<F>)
StantonIntermediatesUCCSD<F>::getRabij(){

  calculateIntermediates();

  if (Rabij) return Rabij;

  //fbd136af1cf7b395e33a6a2040e36b92d7a18e27  -
  Rabij = NEW(CTF::Tensor<F>, *Vabij);

  // P(ab) * Taeij ( Fbe - 0.5 Tbm Fme)
  (*Rabij)["abij"] += (1.0) * (*Tabij)["aeij"] * (*Fae)["be"];
  (*Rabij)["abij"] += (- 1.0) * (*Tabij)["beij"] * (*Fae)["ae"];
  (*Rabij)["abij"] += (- 0.5) * (*Tabij)["aeij"] * (*Tai)["bm"] * (*Fme)["me"];
  (*Rabij)["abij"] += (+ 0.5) * (*Tabij)["beij"] * (*Tai)["am"] * (*Fme)["me"];

  // P(ij) * Tabim ( Fmj + 0.5 Tej Fme)
  (*Rabij)["abij"] += (- 1.0) * (*Tabij)["abim"] * (*Fmi)["mj"];
  (*Rabij)["abij"] += (+ 1.0) * (*Tabij)["abjm"] * (*Fmi)["mi"];
  (*Rabij)["abij"] += (- 0.5) * (*Tabij)["abim"] * (*Tai)["ej"] * (*Fme)["me"];
  (*Rabij)["abij"] += (+ 0.5) * (*Tabij)["abjm"] * (*Tai)["ei"] * (*Fme)["me"];

  (*Rabij)["abij"] += (0.5) * (*Tau_abij)["abmn"] * (*Wijkl)["mnij"];

  //dda7881e81095a6a21011833e88cc972ff890456  -
  //102d94b283918c672a8c9ce73ef689f2f8740661  -
  (*Rabij)["abij"] += (0.5) * (*Tau_abij)["efij"] * (*Wabcd)["abef"];

  //ec0d590a53887b24c495ca90df46ee5782c62515  -
  //e7b0e03e493d603550574c6f4015d16f994eb184  -
  //cff259570ee87e7b824936881b03dc0dc3e55c80  -
  //bcaaa0ebd5951466bf1b3592a990c241c26c7b4e  -
  //aa552987b73a1b0556efa1ab1bb7edf26d5ed58d  -
  // P-ij * P-ab
  (*Rabij)["abij"] += (  1.0) * (*Tabij)["aeim"] * (*Wiabj)["mbej"];
  // -Pij
  (*Rabij)["abij"] += (- 1.0) * (*Tabij)["aejm"] * (*Wiabj)["mbei"];
  // -Pab
  (*Rabij)["abij"] += (- 1.0) * (*Tabij)["beim"] * (*Wiabj)["maej"];
  //  Pij * Pab
  (*Rabij)["abij"] += (  1.0) * (*Tabij)["bejm"] * (*Wiabj)["maei"];

  //896383ad3db1c77dc734c836141429c67784e118  -
  // P-ij * P-ab
  (*Rabij)["abij"] += (- 1.0) * (*Tai)["ei"] * (*Tai)["am"] * (*Viabj)["mbej"];
  // +Pij
  (*Rabij)["abij"] += (+ 1.0) * (*Tai)["ej"] * (*Tai)["am"] * (*Viabj)["mbei"];
  // +Pab
  (*Rabij)["abij"] += (+ 1.0) * (*Tai)["ei"] * (*Tai)["bm"] * (*Viabj)["maej"];
  //  - Pij * Pab
  (*Rabij)["abij"] += (- 1.0) * (*Tai)["ej"] * (*Tai)["bm"] * (*Viabj)["maei"];

  //e035d48a19d337a004e70549c1a490cba713f419  -
  (*Rabij)["abij"] += (  1.0) * (*Tai)["ei"] * (*Vabci)["abej"];
  // - Pij
  (*Rabij)["abij"] += (- 1.0) * (*Tai)["ej"] * (*Vabci)["abei"];

  //b7a00eb9acb88bcf1e3fc73992e621beb09b39f2  -
  (*Rabij)["abij"] += (- 1.0) * (*Tai)["am"] * (*Viajk)["mbij"];
  // + Pab
  (*Rabij)["abij"] += (+ 1.0) * (*Tai)["bm"] * (*Viajk)["maij"];

  return Rabij;

}

// instantiate
template class StantonIntermediatesUCCSD<cc4s::complex>;
template class StantonIntermediatesUCCSD<double>;
