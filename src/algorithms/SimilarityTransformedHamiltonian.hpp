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
      CTF::Tensor<F> *Vabij_ = NULL,
      bool withIntermediates_ = true,
      Dressing dressing_ = Dressing(CCSD)
    );
    virtual ~SimilarityTransformedHamiltonian();

    // dressing tensor setters
    void setTai(CTF::Tensor<F> *Tai_) { Tai = Tai_; }
    void setTabij(CTF::Tensor<F> *Tabij_) { Tabij = Tabij_; }
    void setTabcijk(CTF::Tensor<F> *Tabcijk_) { Tabcijk = Tabcijk_; }

    // ccsd fok vectors
    SDFockVector<F> rightApplyIntermediates(SDFockVector<F> &v);
    SDFockVector<F> rightApplyHirata(SDFockVector<F> &v);
    SDFockVector<F> rightApply(SDFockVector<F> &v);
    SDFockVector<F> leftApplyHirata(SDFockVector<F> &v);

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

    // three body
    PTR(CTF::Tensor<F>) getABCIJK();

    PTR(CTF::Tensor<F>) getTauABIJ();

    void useStantonIntermediatesUCCSD(bool s) {
      _useStantonIntermediatesUCCSD = s;
    }
    bool useStantonIntermediatesUCCSD() {
      return _useStantonIntermediatesUCCSD;
    }

  private:

    // one body
    PTR(CTF::Tensor<F>)  Wij, Wab, Wia, Wai;
    // two body
    PTR(CTF::Tensor<F>)  Wabij, Wijab, Wabcd, Wabci, Waibc,
                         Wiabj, Wiajk, Wijka, Wijkl;
    // three body
    PTR(CTF::Tensor<F>) Wabcijk;

    PTR(CTF::Tensor<F>)  Tau_abij;

    CTF::Tensor<F> *Tai, *Tabij, *Tabcijk;
    CTF::Tensor<F> *Fij, *Fab, *Fia;
    CTF::Tensor<F> *Vabcd, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk,
                   *Vabic, *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci,
                   *Vabij;

    bool _useStantonIntermediatesUCCSD = false;
    PTR(StantonIntermediatesUCCSD<F>) stantonIntermediatesUccsd;
    PTR(StantonIntermediatesUCCSD<F>) getStantonIntermediatesUCCSD();

    bool withIntermediates;

    Dressing dressing;

    int No, Nv;

  };

}

#endif

