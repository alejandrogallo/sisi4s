#ifndef EOM_1RDM_DENSITY_MATRIX_DEFINED
#define EOM_1RDM_DENSITY_MATRIX_DEFINED

#include <iostream>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>
#include <util/Log.hpp>
#include <math/FockVector.hpp>

namespace cc4s {

  template<typename F>
  class OneBodyReducedDensityMatrix {
    public:

    ~OneBodyReducedDensityMatrix() {}
    virtual PTR(CTF::Tensor<F>) getIA() = 0;
    virtual PTR(CTF::Tensor<F>) getAI() = 0;
    virtual PTR(CTF::Tensor<F>) getIJ() = 0;
    virtual PTR(CTF::Tensor<F>) getAB() = 0;

  };

  template<typename F>
  class EomOneBodyReducedDensityMatrix: public OneBodyReducedDensityMatrix<F> {

  public:


    EomOneBodyReducedDensityMatrix(
      CTF::Tensor<F> *Tai_, CTF::Tensor<F> *Tabij_,
      const FockVector<F> *L_, const FockVector<F> *R_
    ): Tai(Tai_), Tabij(Tabij_), L(L_), R(R_){
      LOG(0, "OneBodyRDM") << "Calculating eom ccsd 1rdm" << std::endl;
    }

    EomOneBodyReducedDensityMatrix(
      CTF::Tensor<F> *Tai_, CTF::Tensor<F> *Tabij_, CTF::Tensor<F> *Tabcijk_,
      const FockVector<F> *L_, const FockVector<F> *R_
    ): Tai(Tai_), Tabij(Tabij_), Tabcijk(Tabcijk_), L(L_), R(R_){
      LOG(0, "OneBodyRDM") << "Calculating eom ccsdt 1rdm" << std::endl;
    }

    PTR(CTF::Tensor<F>) getAI() {
      if (Rai) { return Rai; }
      Rai = NEW(CTF::Tensor<F>, *Tai);

      (*Rai)["ai"]  = 0;
      (*Rai)["ai"] += (*L->get(0))["ke"] * (*R->get(1))["eaki"];
      (*Rai)["ai"] += (*L->get(0))["ke"] * (*R->get(0))["ak"] * (*Tai)["ei"];
      (*Rai)["ai"] += (*L->get(0))["ke"] * (*R->get(0))["ei"] * (*Tai)["ak"];
      (*Rai)["ai"] += (-0.5) * (*L->get(1))["kled"] * (*R->get(0))["di"] * (*Tabij)["eakl"];
      (*Rai)["ai"] += (-0.5) * (*L->get(1))["kled"] * (*R->get(0))["al"] * (*Tabij)["edki"];
      (*Rai)["ai"] += (-0.5) * (*L->get(1))["kled"] * (*Tai)["di"] * (*R->get(1))["eakl"];
      (*Rai)["ai"] += (-0.5) * (*L->get(1))["kled"] * (*Tai)["al"] * (*R->get(1))["edki"];

      return Rai;
    }

    PTR(CTF::Tensor<F>) getIJ() {
      if (Rij) { return Rij; }
      const int No(Tai->lens[1]);
      const int oo[] = {No, No};
      const int syms[] = {NS, NS};
      Rij = NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "Rij");

      (*Rij)["ij"]  = 0;
      (*Rij)["ij"] += (*L->get(0))["je"] * (*R->get(0))["ei"];
      (*Rij)["ij"] += 0.5 * (*L->get(1))["kjed"] * (*R->get(1))["edki"];
      (*Rij)["ij"] += (*L->get(1))["kjed"] * (*R->get(0))["ek"] * (*Tai)["di"];
      // This is not in the paper
      (*Rij)["ij"] += (*L->get(1))["kjed"] * (*R->get(0))["di"] * (*Tai)["ek"];

      return Rij;
    }

    PTR(CTF::Tensor<F>) getAB() {
      if (Rab) { return Rab; }
      const int Nv(Tai->lens[0]);
      const int vv[] = {Nv, Nv};
      const int syms[] = {NS, NS};
      Rab = NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "Rab");

      (*Rab)["ab"]  = 0;
      (*Rab)["ab"] += (-1.0) * (*L->get(0))["ka"] * (*R->get(0))["bk"];
      (*Rab)["ab"] += (-0.5) * (*L->get(1))["klea"] * (*R->get(1))["ebkl"];
      (*Rab)["ab"] += (-1.0) * (*L->get(1))["klea"] * (*R->get(0))["ek"] * (*Tai)["bl"];
      // This is not in the paper
      (*Rab)["ab"] += (-1.0) * (*L->get(1))["klea"] * (*R->get(0))["bl"] * (*Tai)["ek"];

      return Rab;

    }

    PTR(CTF::Tensor<F>) getIA() {
      if (Ria) { return Ria; }
      const int Nv(Tai->lens[0]), No(Tai->lens[1]);
      const int ov[] = {No, Nv};
      const int syms[] = {NS, NS};

      Ria = NEW(CTF::Tensor<F>, 2, ov, syms, *Cc4s::world, "Ria");

      (*Ria)["ia"]  = 0;
      // this is 0 because r0 is 0
      //(*Ria)["ia"] += (*L->get(0))["ia"];
      (*Ria)["ia"] += (*L->get(1))["oifa"] * (*R->get(0))["fo"];
      // this is 0 because r0 is 0
      //(*Ria)["ia"] += (*L->get(1))["oifa"] * Tai["fo"];

      return Ria;
    }

  private:

    PTR(CTF::Tensor<F>) Rai, Ria, Rij, Rab;
    CTF::Tensor<F> *Tai, *Tabij, *Tabcijk;
    const FockVector<F> *L, *R;

  };

}

#endif