#ifndef SIM_TRANS_HAMILTONIAN_DEFINED
#define SIM_TRANS_HAMILTONIAN_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/StantonIntermediatesUCCSD.hpp>
#include <math/FockVector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>

namespace cc4s {

  template <typename F=complex>
  class SimilarityTransformedHamiltonian {
  public:
    typedef SimilarityTransformedHamiltonian STH;

    /*! \enum Dressing
     *
     *  The kind of dressing of the Similarity transformed hamiltonian,
     *  i.e., the nature of T, if T is the coupled cluster excitation
     *  operator, then some matrices elements of H will be identically 0.
     */
    enum Dressing {
      CCSD,
      CCSDT,
      NONE, // Some terms will be then set to zero instead of multiplying by 0
      GENERAL, // For a general T not fulfiling any particular criterium
    };

    SimilarityTransformedHamiltonian(
      int No, int Nv,
      bool withIntermediates = true
    );
    virtual ~SimilarityTransformedHamiltonian();

    // ccsd fok vectors
    SDFockVector<F> rightApplyIntermediates(SDFockVector<F> &v);
    SDFockVector<F> rightApplyHirata(SDFockVector<F> &v);
    SDFockVector<F> rightApply(SDFockVector<F> &v);
    SDFockVector<F> leftApplyHirata(SDFockVector<F> &v);
    SDFockVector<F> leftApplyIntermediates(SDFockVector<F> &v);
    SDFockVector<F> leftApply(SDFockVector<F> &v);

    // ccsdt fok vectors
    SDTFockVector<F> rightApplyHirata(SDTFockVector<F> &v);
    SDTFockVector<F> rightApply(SDTFockVector<F> &v);

    // One body
    PTR(CTF::Tensor<F>) getIJ();
    PTR(CTF::Tensor<F>) getAB();
    PTR(CTF::Tensor<F>) getAI();
    PTR(CTF::Tensor<F>) getIA();

    // Two body
    PTR(CTF::Tensor<F>) getABIJ();
    PTR(CTF::Tensor<F>) getIJAB();
    PTR(CTF::Tensor<F>) getABCD();
    PTR(CTF::Tensor<F>) getABCI();
    PTR(CTF::Tensor<F>) getAIBC();
    PTR(CTF::Tensor<F>) getIABJ();
    PTR(CTF::Tensor<F>) getIAJK();
    PTR(CTF::Tensor<F>) getIJKA();
    PTR(CTF::Tensor<F>) getIJKL();

    // dressing tensor setters
    STH& setTai(CTF::Tensor<F> *t) { Tai = t; return *this;}
    STH& setTabij(CTF::Tensor<F> *t) { Tabij = t; return *this;}
    STH& setTabcijk(CTF::Tensor<F> *t) { Tabcijk = t; return *this;}

    // V amplitudes setters
    STH& setFij(CTF::Tensor<F> *t) { Fij = t; return *this; }
    STH& setFab(CTF::Tensor<F> *t) { Fab = t; return *this; }
    STH& setFia(CTF::Tensor<F> *t) { Fia = t; return *this; }
    STH& setVabcd(CTF::Tensor<F> *t) { Vabcd = t; return *this; }
    STH& setViajb(CTF::Tensor<F> *t) { Viajb = t; return *this; }
    STH& setVijab(CTF::Tensor<F> *t) { Vijab = t; return *this; }
    STH& setVijkl(CTF::Tensor<F> *t) { Vijkl = t; return *this; }
    STH& setVijka(CTF::Tensor<F> *t) { Vijka = t; return *this; }
    STH& setViabc(CTF::Tensor<F> *t) { Viabc = t; return *this; }
    STH& setViajk(CTF::Tensor<F> *t) { Viajk = t; return *this; }
    STH& setVabic(CTF::Tensor<F> *t) { Vabic = t; return *this; }
    STH& setVaibc(CTF::Tensor<F> *t) { Vaibc = t; return *this; }
    STH& setVaibj(CTF::Tensor<F> *t) { Vaibj = t; return *this; }
    STH& setViabj(CTF::Tensor<F> *t) { Viabj = t; return *this; }
    STH& setVijak(CTF::Tensor<F> *t) { Vijak = t; return *this; }
    STH& setVaijb(CTF::Tensor<F> *t) { Vaijb = t; return *this; }
    STH& setVabci(CTF::Tensor<F> *t) { Vabci = t; return *this; }
    STH& setVabij(CTF::Tensor<F> *t) { Vabij = t; return *this; }

    STH& setWithIntermediates(bool t) {withIntermediates = t; return *this;}
    STH& setDressing(Dressing d) {dressing = d; return *this;}

    // three body
    PTR(CTF::Tensor<F>) getABCIJK();

    PTR(CTF::Tensor<F>) getTauABIJ();

    STH& useStantonIntermediatesUCCSD(bool s) {
      _useStantonIntermediatesUCCSD = s;
      return *this;
    }
    bool useStantonIntermediatesUCCSD() {
      return _useStantonIntermediatesUCCSD;
    }

  private:

    bool _useStantonIntermediatesUCCSD = false;
    PTR(StantonIntermediatesUCCSD<F>) stantonIntermediatesUccsd;
    PTR(StantonIntermediatesUCCSD<F>) getStantonIntermediatesUCCSD();

    int No, Nv;
    bool withIntermediates;
    Dressing dressing;

    //
    // Resources that should be destroyed after the class gets destroyed
    //
    // one body
    PTR(CTF::Tensor<F>)  Wij, Wab, Wia, Wai;
    // two body
    PTR(CTF::Tensor<F>)  Wabij, Wijab, Wabcd, Wabci, Waibc,
                         Wiabj, Wiajk, Wijka, Wijkl;
    // three body
    PTR(CTF::Tensor<F>) Wabcijk;

    PTR(CTF::Tensor<F>)  Tau_abij;


    //
    // External resources that should not be cleaned up after
    // Hamiltonian class gets destroyed
    //
    CTF::Tensor<F> *Tai, *Tabij, *Tabcijk;
    CTF::Tensor<F> *Fij, *Fab, *Fia;
    CTF::Tensor<F> *Vabcd, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk,
                   *Vabic, *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci,
                   *Vabij;

  };

}

#endif
