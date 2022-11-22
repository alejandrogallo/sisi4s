#ifndef SIM_TRANS_HAMILTONIAN_DEFINED
#define SIM_TRANS_HAMILTONIAN_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/StantonIntermediatesUCCSD.hpp>
#include <math/FockVector.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>

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

  SimilarityTransformedHamiltonian(int No_, int Nv_)
      : No(No_)
      , Nv(Nv_){};
  ~SimilarityTransformedHamiltonian(){};

  // RPA, singles
  SFockVector<F> rightApplyHirata_RPA(SFockVector<F> &v);

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

  // ip and ea ccsd
  SDFockVector<F> rightApply_CCSD_IP(SDFockVector<F> &);
  SDFockVector<F> rightApplyHirata_CCSD_IP(SDFockVector<F> &);
  SDFockVector<F> rightApplyIntermediates_CCSD_IP(SDFockVector<F> &);

  SDFockVector<F> rightApply_CCSD_EA(SDFockVector<F> &);
  SDFockVector<F> rightApplyHirata_CCSD_EA(SDFockVector<F> &);
  SDFockVector<F> rightApplyIntermediates_CCSD_EA(SDFockVector<F> &);

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

  // dressing tensor setters
  STH &setTai(Tensor<F> *t) {
    Tai = t;
    return *this;
  }
  STH &setTabij(Tensor<F> *t) {
    Tabij = t;
    return *this;
  }
  STH &setTabcijk(Tensor<F> *t) {
    Tabcijk = t;
    return *this;
  }
  STH &setTabcdijkl(Tensor<F> *t) {
    Tabcdijkl = t;
    return *this;
  }

  // V amplitudes setters
  STH &setFij(Tensor<F> *t) {
    Fij = t;
    return *this;
  }
  STH &setFab(Tensor<F> *t) {
    Fab = t;
    return *this;
  }
  STH &setFia(Tensor<F> *t) {
    Fia = t;
    return *this;
  }
  STH &setVabcd(Tensor<F> *t) {
    Vabcd = t;
    return *this;
  }
  STH &setViajb(Tensor<F> *t) {
    Viajb = t;
    return *this;
  }
  STH &setVijab(Tensor<F> *t) {
    Vijab = t;
    return *this;
  }
  STH &setVijkl(Tensor<F> *t) {
    Vijkl = t;
    return *this;
  }
  STH &setVijka(Tensor<F> *t) {
    Vijka = t;
    return *this;
  }
  STH &setViabc(Tensor<F> *t) {
    Viabc = t;
    return *this;
  }
  STH &setViajk(Tensor<F> *t) {
    Viajk = t;
    return *this;
  }
  STH &setVabic(Tensor<F> *t) {
    Vabic = t;
    return *this;
  }
  STH &setVaibc(Tensor<F> *t) {
    Vaibc = t;
    return *this;
  }
  STH &setVaibj(Tensor<F> *t) {
    Vaibj = t;
    return *this;
  }
  STH &setViabj(Tensor<F> *t) {
    Viabj = t;
    return *this;
  }
  STH &setVijak(Tensor<F> *t) {
    Vijak = t;
    return *this;
  }
  STH &setVaijb(Tensor<F> *t) {
    Vaijb = t;
    return *this;
  }
  STH &setVabci(Tensor<F> *t) {
    Vabci = t;
    return *this;
  }
  STH &setVabij(Tensor<F> *t) {
    Vabij = t;
    return *this;
  }

  // coulomb bertex setter
  STH &setGammaGqr(Tensor<sisi4s::complex> *t) {
    GammaGqr = t;
    return *this;
  }

  STH &setRightApplyIntermediates(bool t) {
    useRightApplyIntermediates = t;
    return *this;
  }
  STH &setDressing(Dressing d) {
    dressing = d;
    return *this;
  }

  PTR(Tensor<F>) getTauABIJ();

  STH &useStantonIntermediatesUCCSD(bool s) {
    _useStantonIntermediatesUCCSD = s;
    return *this;
  }
  bool useStantonIntermediatesUCCSD() { return _useStantonIntermediatesUCCSD; }

  std::string getAbbreviation() const { return "STH"; }

private:
  bool _useStantonIntermediatesUCCSD = false;
  PTR(StantonIntermediatesUCCSD<F>) stantonIntermediatesUccsd;
  PTR(StantonIntermediatesUCCSD<F>) getStantonIntermediatesUCCSD();

  int No, Nv;
  bool useRightApplyIntermediates;
  Dressing dressing;

  //
  // Resources that should be destroyed after the class gets destroyed
  //
  PTR(Tensor<F>)
  // one body
  Wij, Wab, Wia,
      Wai

      // two body
      ,
      Wabij, Wijab, Wabcd, Wabci, Waibc, Wiabj, Wiajk, Wijka,
      Wijkl

      // three body
      ,
      Wabcijk

      // intermediate quantities
      ,
      Tau_abij;

  //
  // External resources that should not be cleaned up after
  // Hamiltonian class gets destroyed
  //
  Tensor<F>
      // T amplitudes
      *Tai = nullptr,
      *Tabij = nullptr, *Tabcijk = nullptr,
      *Tabcdijkl = nullptr

      // Fock matrices
      ,
      *Fij, *Fab,
      *Fia = nullptr

      // Coulomb integrals
      ,
      *Vabcd = nullptr, *Viajb, *Vijab, *Vijkl, *Vijka, *Viabc, *Viajk, *Vabic,
      *Vaibc, *Vaibj, *Viabj, *Vijak, *Vaijb, *Vabci, *Vabij;

  // coulomb vertex
  Tensor<sisi4s::complex> *GammaGqr = nullptr;
};

} // namespace sisi4s

#endif
