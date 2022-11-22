#ifndef STANTON_UCCSD_INTERMEDIATES_DEFINED
#define STANTON_UCCSD_INTERMEDIATES_DEFINED

#include <util/SharedPointer.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>

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
template <typename F = complex>
class StantonIntermediatesUCCSD {
public:
  /**
   * \brief Get the singles part of the residuum
   */
  PTR(Tensor<F>) getRai();
  /**
   * \brief Get the doubles part of the residuum
   */
  PTR(Tensor<F>) getRabij();

  StantonIntermediatesUCCSD() {}

  // tensor setters, it's more convenient to set them like so
  // than in a constructor, since there are so many tensors
  void setTai(Tensor<F> *t) { Tai = t; }
  void setTabij(Tensor<F> *t) { Tabij = t; }
  void setVabcd(Tensor<F> *t) { Vabcd = t; }
  void setViajb(Tensor<F> *t) { Viajb = t; }
  void setVijab(Tensor<F> *t) { Vijab = t; }
  void setVijkl(Tensor<F> *t) { Vijkl = t; }
  void setVijka(Tensor<F> *t) { Vijka = t; }
  void setViabc(Tensor<F> *t) { Viabc = t; }
  void setViajk(Tensor<F> *t) { Viajk = t; }
  void setVabic(Tensor<F> *t) { Vabic = t; }
  void setVaibc(Tensor<F> *t) { Vaibc = t; }
  void setVaibj(Tensor<F> *t) { Vaibj = t; }
  void setViabj(Tensor<F> *t) { Viabj = t; }
  void setVijak(Tensor<F> *t) { Vijak = t; }
  void setVaijb(Tensor<F> *t) { Vaijb = t; }
  void setVabci(Tensor<F> *t) { Vabci = t; }
  void setVabij(Tensor<F> *t) { Vabij = t; }
  void setFij(Tensor<F> *t) { Fij = t; }
  void setFab(Tensor<F> *t) { Fab = t; }
  void setFia(Tensor<F> *t) { Fia = t; }

private:
  void calculateIntermediates();
  void calculateOneBodyIntermediates();
  void checkInputs();

  // Intermediate quantities
  PTR(Tensor<F>)
  Fae, Fmi, Fme, Tau_abij, TildeTau_abij, Wijkl, Wabcd, Wiabj;

  // Output quantities
  PTR(Tensor<F>) Rai, Rabij;

  // Input quantities
  Tensor<F> *Tai, *Tabij;
  Tensor<F> *Fij, *Fab, *Fia = nullptr;
  Tensor<F> *Vabcd, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk, *Vabic,
      *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci, *Vabij;
};

} // namespace sisi4s

#endif
