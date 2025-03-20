#ifndef SIM_TRANS_HAMILTONIAN_DEFINED
#define SIM_TRANS_HAMILTONIAN_DEFINED

#include <equations/CoulombIntegrals.hpp>
#include <equations/StantonIntermediatesUCCSD.hpp>
#include <math/FockVector.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>

// Utility Macros

#define _DEFINE_SETTER(type, name, default)                                    \
  SimilarityTransformedHamiltonian &set##name(type t) {                        \
    name = t;                                                                  \
    return *this;                                                              \
  }                                                                            \
  type name = default

#define _DEFINE_V_SETTER(type, name, indices)                                  \
  SimilarityTransformedHamiltonian &set##name(type t) {                        \
    if (Vpqrs == nullptr) { Vpqrs = new CoulombIntegrals<F>(); }               \
    Vpqrs->with_##indices(t);                                                  \
    return *this;                                                              \
  }

#define _MAKE_WITH_FUNCTION(type, name, default)                               \
  SimilarityTransformedHamiltonian &with##name(type const &v) {                \
    _with##name = v;                                                           \
    return *this;                                                              \
  }                                                                            \
  type &with##name() { return _with##name; }                                   \
  type _with##name = default

namespace sisi4s {

template <typename F = complex>
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
    RPA,
    GENERAL, // For a general T not fulfiling any particular criterium
  };

  SimilarityTransformedHamiltonian(const int No_, const int Nv_)
      : No(No_)
      , Nv(Nv_){};
  SimilarityTransformedHamiltonian(const int No_,
                                   const int Nv_,
                                   CoulombIntegrals<F> *v)
      : No(No_)
      , Nv(Nv_)
      , Vpqrs(v){};
  ~SimilarityTransformedHamiltonian(){};

  // RPA, singles
  SFockVector<F> right_apply_hirata_RPA(SFockVector<F> &v);

  // CISD
  FockVector<F> right_apply_CISD(FockVector<F> &R);

  // ccsd fok vectors
  SDFockVector<F> right_apply_Intermediates(SDFockVector<F> &v);
  SDFockVector<F> right_apply_hirata(SDFockVector<F> &v);
  SDFockVector<F> right_apply(SDFockVector<F> &v);
  SDFockVector<F> leftApply_hirata(SDFockVector<F> &v);
  SDFockVector<F> leftApplyIntermediates(SDFockVector<F> &v);
  SDFockVector<F> leftApply(SDFockVector<F> &v);

  // ccsdt fok vectors
  SDTFockVector<F> right_apply_hirata(SDTFockVector<F> &v);
  SDTFockVector<F> right_apply(SDTFockVector<F> &v);

  // ip and ea ccsd
  SDFockVector<F> right_apply_CCSD_IP(SDFockVector<F> &);
  SDFockVector<F> right_apply_hirata_CCSD_IP(SDFockVector<F> &);
  SDFockVector<F> right_apply_Intermediates_CCSD_IP(SDFockVector<F> &);

  SDFockVector<F> right_apply_CCSD_EA(SDFockVector<F> &);
  SDFockVector<F> right_apply_hirata_CCSD_EA(SDFockVector<F> &);
  SDFockVector<F> right_apply_Intermediates_CCSD_EA(SDFockVector<F> &);

  // block2
  FockVector<F> right_apply_CCSDT_IP_BLOCK2(FockVector<F> &);
  FockVector<F> right_apply_CCSDT_EA_BLOCK2(FockVector<F> &);

  // Structure factor
  struct StructureFactor {
    double energy;
    Tensor<F> S;
  };
  struct StructureFactorSettings {
    bool onlySingles;
    bool onlyDoubles;
    bool hartreeInOneBody;
    bool fockInOneBody;
  };
  StructureFactor structureFactor(SDFockVector<F> &t,
                                  const StructureFactorSettings &s = {
                                      .onlySingles = false,
                                      .onlyDoubles = false,
                                      .hartreeInOneBody = false,
                                      .fockInOneBody = false});

  //
  // External resources that should not be cleaned up after
  // Hamiltonian class gets destroyed
  //
  // One body
  PTR(Tensor<F>) getIJ();
  PTR(Tensor<F>) getAB();
  PTR(Tensor<F>) getAI();
  PTR(Tensor<F>) getAI_RPA();
  PTR(Tensor<F>) getIA();

  // Two body
  PTR(Tensor<F>) getABIJ();
  PTR(Tensor<F>) getABIJ_RPA();
  PTR(Tensor<F>) getIJAB();
  PTR(Tensor<F>) getABCD();
  PTR(Tensor<F>) getABCI();
  PTR(Tensor<F>) getAIBC();
  PTR(Tensor<F>) getIABJ();
  PTR(Tensor<F>) getIAJK();
  PTR(Tensor<F>) getIJKA();
  PTR(Tensor<F>) getIJKL();

  // three body
  PTR(Tensor<F>) getABCIJK();

  // four body
  PTR(Tensor<F>) getABCDIJKL();

  // t amplitudes
  _DEFINE_SETTER(Tensor<F> *, Tai, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Tabij, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Tabcijk, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Tabcdijkl, nullptr);

  // one body part
  _DEFINE_SETTER(Tensor<F> *, Fij, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Fab, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Fia, nullptr);

  // Coulomb Integrals
  _DEFINE_V_SETTER(Tensor<F> *, Vabcd, pppp);
  _DEFINE_V_SETTER(Tensor<F> *, Viajb, hphp);
  _DEFINE_V_SETTER(Tensor<F> *, Vijab, hhpp);
  _DEFINE_V_SETTER(Tensor<F> *, Vijkl, hhhh);
  _DEFINE_V_SETTER(Tensor<F> *, Vijka, hhhp);
  _DEFINE_V_SETTER(Tensor<F> *, Viabc, hppp);
  _DEFINE_V_SETTER(Tensor<F> *, Viajk, hphh);
  _DEFINE_V_SETTER(Tensor<F> *, Vabic, pphp);
  _DEFINE_V_SETTER(Tensor<F> *, Vaibc, phpp);
  _DEFINE_V_SETTER(Tensor<F> *, Vaibj, phph);
  _DEFINE_V_SETTER(Tensor<F> *, Viabj, hpph);
  _DEFINE_V_SETTER(Tensor<F> *, Vijak, hhph);
  _DEFINE_V_SETTER(Tensor<F> *, Vaijk, phhh);
  _DEFINE_V_SETTER(Tensor<F> *, Vaijb, phhp);
  _DEFINE_V_SETTER(Tensor<F> *, Vabci, ppph);
  _DEFINE_V_SETTER(Tensor<F> *, Vabij, pphh);

  _DEFINE_SETTER(Tensor<F> *, VVaijb, nullptr);
  _DEFINE_SETTER(Tensor<F> *, VViabc, nullptr);
  _DEFINE_SETTER(Tensor<F> *, VVijka, nullptr);
  _DEFINE_SETTER(Tensor<F> *, VVijab, nullptr);

  // coulomb vertex setter
  _DEFINE_SETTER(Tensor<sisi4s::complex> *, GammaGqr, nullptr);

  STH &setDressing(Dressing d) {
    dressing = d;
    return *this;
  }

  PTR(Tensor<F>) getTauABIJ();

  std::string getAbbreviation() const { return "STH"; }

  _MAKE_WITH_FUNCTION(bool, _right_apply_intermediates, false);
  _MAKE_WITH_FUNCTION(bool, StantonIntermediatesUCCSD, false);
  _MAKE_WITH_FUNCTION(bool, RingCCSDT, false);
  _MAKE_WITH_FUNCTION(bool, CISD, false);
  _MAKE_WITH_FUNCTION(bool, CIS, false);

private:
  PTR(StantonIntermediatesUCCSD<F>) stantonIntermediatesUccsd;
  PTR(StantonIntermediatesUCCSD<F>) getStantonIntermediatesUCCSD();

  int No, Nv;
  Dressing dressing;

  CoulombIntegrals<F> *Vpqrs = nullptr;

  //
  // Resources that should be destroyed after the class gets destroyed
  //
  PTR(Tensor<F>)
  // one body
  Wij, Wab, Wia, Wai,

      // two body
      Wabij, Wijab, Wabcd, Wabci, Waibc, Wiabj, Wiajk, Wijka, Wijkl,

      // three body
      Wabcijk,

      // four body
      Wabcdijkl,

      // intermediate quantities
      Tau_abij;
};

} // namespace sisi4s

#undef _DEFINE_SETTER
#undef _DEFINE_V_SETTER
#undef _MAKE_WITH_FUNCTION

#endif
