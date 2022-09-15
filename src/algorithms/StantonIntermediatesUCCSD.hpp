#ifndef STANTON_UCCSD_INTERMEDIATES_DEFINED
#define STANTON_UCCSD_INTERMEDIATES_DEFINED

#include <util/SharedPointer.hpp>
#include <math/Complex.hpp>
#include <util/CTF.hpp>

namespace sisi4s {

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
    void setTai(CTF::Tensor<F> *t) { Tai = t; }
    void setTabij(CTF::Tensor<F> *t) { Tabij = t; }
    void setVabcd(CTF::Tensor<F> *t) { Vabcd = t; }
    void setViajb(CTF::Tensor<F> *t) { Viajb = t; }
    void setVijab(CTF::Tensor<F> *t) { Vijab = t; }
    void setVijkl(CTF::Tensor<F> *t) { Vijkl = t; }
    void setVijka(CTF::Tensor<F> *t) { Vijka = t; }
    void setViabc(CTF::Tensor<F> *t) { Viabc = t; }
    void setViajk(CTF::Tensor<F> *t) { Viajk = t; }
    void setVabic(CTF::Tensor<F> *t) { Vabic = t; }
    void setVaibc(CTF::Tensor<F> *t) { Vaibc = t; }
    void setVaibj(CTF::Tensor<F> *t) { Vaibj = t; }
    void setViabj(CTF::Tensor<F> *t) { Viabj = t; }
    void setVijak(CTF::Tensor<F> *t) { Vijak = t; }
    void setVaijb(CTF::Tensor<F> *t) { Vaijb = t; }
    void setVabci(CTF::Tensor<F> *t) { Vabci = t; }
    void setVabij(CTF::Tensor<F> *t) { Vabij = t; }
    void setFij(CTF::Tensor<F> *t) { Fij = t; }
    void setFab(CTF::Tensor<F> *t) { Fab = t; }
    void setFia(CTF::Tensor<F> *t) { Fia = t; }

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
