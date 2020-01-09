#ifndef STANTON_UCCSD_INTERMEDIATES_DEFINED
#define STANTON_UCCSD_INTERMEDIATES_DEFINED

#include <util/SharedPointer.hpp>
#include <math/Complex.hpp>
#include <ctf.hpp>

namespace cc4s {

/*
 * This is the code version of:
 * ----------------------------
 * John F. Stanton, Jürgen Gauss, John D. Watts, Rodney J. Bartlett.
 * A direct product decomposition approach for symmetry exploitation in
 * many‐body methods. I. Energy calculations.
 * The Journal of Chemical Physics  1991 doi:10.1063/1.460620
 *
 */
  template <typename F=complex>
  class StantonIntermediatesUCCSD {
  public:
    /**
     * \brief Get the singles part of the residuum
     */
    PTR(CTF::Tensor<F>) getRai();
    /**
     * \brief Get the doubles part of the residuum
     */
    PTR(CTF::Tensor<F>) getRabij();

    StantonIntermediatesUCCSD() {}

    // tensor setters, it's more convenient to set them like so
    // than in a constructor, since there are so many tensors
    void setTai(CTF::Tensor<F> *Tai_) { Tai = Tai_; }
    void setTabij(CTF::Tensor<F> *Tabij_) { Tabij = Tabij_; }
    void setVabcd(CTF::Tensor<F> *Vabcd_) { Vabcd = Vabcd_; }
    void setViajb(CTF::Tensor<F> *Viajb_) { Viajb = Viajb_; }
    void setVijab(CTF::Tensor<F> *Vijab_) { Vijab = Vijab_; }
    void setVijkl(CTF::Tensor<F> *Vijkl_) { Vijkl = Vijkl_; }
    void setVijka(CTF::Tensor<F> *Vijka_) { Vijka = Vijka_; }
    void setViabc(CTF::Tensor<F> *Viabc_) { Viabc = Viabc_; }
    void setViajk(CTF::Tensor<F> *Viajk_) { Viajk = Viajk_; }
    void setVabic(CTF::Tensor<F> *Vabic_) { Vabic = Vabic_; }
    void setVaibc(CTF::Tensor<F> *Vaibc_) { Vaibc = Vaibc_; }
    void setVaibj(CTF::Tensor<F> *Vaibj_) { Vaibj = Vaibj_; }
    void setViabj(CTF::Tensor<F> *Viabj_) { Viabj = Viabj_; }
    void setVijak(CTF::Tensor<F> *Vijak_) { Vijak = Vijak_; }
    void setVaijb(CTF::Tensor<F> *Vaijb_) { Vaijb = Vaijb_; }
    void setVabci(CTF::Tensor<F> *Vabci_) { Vabci = Vabci_; }
    void setVabij(CTF::Tensor<F> *Vabij_) { Vabij = Vabij_; }
    void setFij(CTF::Tensor<F> *Fij_) { Fij = Fij_; }
    void setFab(CTF::Tensor<F> *Fab_) { Fab = Fab_; }
    void setFia(CTF::Tensor<F> *Fia_) { Fia = Fia_; }

  private:

    void calculateIntermediates();
    void calculateOneBodyIntermediates();
    void checkInputs();

    // Intermediate quantities
    PTR(CTF::Tensor<F>)
      Fae, Fmi, Fme, Tau_abij, TildeTau_abij, Wijkl, Wabcd, Wiabj;

    // Output quantities
    PTR(CTF::Tensor<F>) Rai, Rabij;

    // Input quantities
    CTF::Tensor<F> *Tai, *Tabij;
    CTF::Tensor<F> *Fij, *Fab, *Fia=nullptr;
    CTF::Tensor<F>
      *Vabcd, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk,
      *Vabic, *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci, *Vabij;
  };

}

#endif
