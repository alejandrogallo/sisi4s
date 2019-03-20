#ifndef CCSD_SIM_TRANS_DEFINED
#define CCSD_SIM_TRANS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>

namespace cc4s {

  template <typename F=complex>
  class CcsdSimilarityTransformedHamiltonian {
  public:
    CcsdSimilarityTransformedHamiltonian(
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
      CTF::Tensor<F> *Vabci_
    );
    virtual ~CcsdSimilarityTransformedHamiltonian();
    virtual void run();

    FockVector<F> rightApplyIntermediates(FockVector<F> &v);
    FockVector<F> rightApplyHirata(FockVector<F> &v);
    FockVector<F> rightApply(FockVector<F> &v);
    FockVector<F> leftApply(FockVector<F> &v);

    // One body
    PTR(CTF::Tensor<F>) getWij();
    CTF::Tensor<F>* getWab();
    CTF::Tensor<F>* getWia();
    CTF::Tensor<F>* getWai();

    // Two body
    CTF::Tensor<F>* getWabij();
    CTF::Tensor<F>* getWijab();
    CTF::Tensor<F>* getWabcd();
    CTF::Tensor<F>* getWabci();
    CTF::Tensor<F>* getWaibc();
    CTF::Tensor<F>* getWiabj();
    CTF::Tensor<F>* getWiajk();
    CTF::Tensor<F>* getWijka();
    CTF::Tensor<F>* getWijkl();

    void
    buildAllIntermediates(bool intermediates);

  private:
    PTR(CTF::Tensor<F>)  Wij, Wab, Wia, Wai;
    PTR(CTF::Tensor<F>)  Wabij, Wijab, Wabcd, Wabci, Waibc,
                         Wiabj, Wiajk, Wijka, Wijkl;

    CTF::Tensor<F> *Tai, *Tabij;
    PTR(CTF::Tensor<F>)  Tau_abij;
    CTF::Tensor<F> *Fij, *Fab, *Fia;
    CTF::Tensor<F> *Vabcd, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk,
                   *Vabic, *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci;

    int No, Nv;

    bool withIntermediates;

  };

  class CcsdSimilarityTransformedHamiltonianFactory: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdSimilarityTransformedHamiltonianFactory);
    CcsdSimilarityTransformedHamiltonianFactory(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdSimilarityTransformedHamiltonianFactory();
    virtual void run();

    template<typename F>
    void run();

  };

}

#endif

