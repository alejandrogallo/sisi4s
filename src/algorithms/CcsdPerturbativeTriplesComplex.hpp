#ifndef CCSD_PERTURBATIVE_TRIPLES_COMPLEX_DEFINED
#define CCSD_PERTURBATIVE_TRIPLES_COMPLEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Permutation.hpp>
#include <util/SlicedCtfTensor.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
class CcsdPerturbativeTriplesComplex : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdPerturbativeTriplesComplex);
  CcsdPerturbativeTriplesComplex(std::vector<Argument> const &argumentList);
  virtual ~CcsdPerturbativeTriplesComplex();
  /**
   * \brief Calculates perturbative triples correction. Routine based on
   * Helgaker book.
   */
  virtual void run();

  /**
   * \brief Dry run for perturbative triples correction based on Helgaker book.
   */
  virtual void dryRun();

private:
  // NOTE: the Dummy template argument is needed to "fully" specialize
  // the inner class CoulombVertex to CoulombVertex<double> or <complex>
  template <typename F, int Dummy = 0>
  class CoulombVertex {};

  template <int Dummy>
  class CoulombVertex<double, Dummy> {
  public:
    CoulombVertex(Tensor<complex> *GammaFab, Tensor<complex> *GammaFai);
    ~CoulombVertex();
    void getDoublesParticleContribution(SlicedCtfTensor<double> &Tabij,
                                        const Map<3> &i,
                                        Tensor<double> &SVabc);

  protected:
    Tensor<double> *realGammaFab, *imagGammaFab;
    SlicedCtfTensor<double> *realGammaFai, *imagGammaFai;
  };

  template <int Dummy>
  class CoulombVertex<complex, Dummy> {
  public:
    CoulombVertex(Tensor<complex> *GammaFab, Tensor<complex> *GammaFai);
    ~CoulombVertex();
    void getDoublesParticleContribution(SlicedCtfTensor<complex> &Tabij,
                                        const Map<3> &i,
                                        Tensor<complex> &DVabc);

  protected:
    Tensor<complex> *conjGammaFab;
    SlicedCtfTensor<complex> *GammaFai;
  };

  template <typename F>
  class Calculator {
  public:
    Calculator(Tensor<F> *Tai,
               Tensor<F> *Tabij,
               Tensor<F> *Vabij,
               Tensor<F> *Valij,
               Tensor<complex> *GammaFab,
               Tensor<complex> *GammaFai,
               Tensor<double> *epsi,
               Tensor<double> *epsa);
    ~Calculator();
    F calculate();
    void addDoublesHoleContribution(const Map<3> &, Tensor<F> &);
    Tensor<F> &getSinglesContribution(const Map<3> &);
    Tensor<F> &getEnergyDenominator(const Map<3> &);

  protected:
    Tensor<F> *SVabc, *DVabc;
    SlicedCtfTensor<F> *Tai, *Tabij, *Tabil;
    SlicedCtfTensor<F> *Vabij, *Valij;
    CoulombVertex<F> Gamma;
    SlicedCtfTensor<F> *epsi;
    Tensor<F> *epsa;
  };
};
} // namespace sisi4s

#endif
