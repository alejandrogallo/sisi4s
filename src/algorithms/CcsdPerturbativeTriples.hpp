#ifndef CCSD_PERTURBATIVE_TRIPLES_DEFINED
#define CCSD_PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Permutation.hpp>
#include <util/SlicedCtfTensor.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
class CcsdPerturbativeTriples : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdPerturbativeTriples);
  CcsdPerturbativeTriples(std::vector<Argument> const &argumentList);
  virtual ~CcsdPerturbativeTriples();
  /**
   * \brief Calculates perturbative triples correction. Routine based on
   * Helgaker book.
   */
  virtual void run();

  /**
   * \brief Dry run for perturbative triples correction based on Helgaker book.
   */
  virtual void dryRun();

protected:
  int No, Nv;
  Tensor<double> *SVabc, *DVabc;
  Tensor<double> *realGammaFab, *imagGammaFab;
  SlicedCtfTensor<> *Tai, *Tabij, *Tabil;
  SlicedCtfTensor<> *Vabij, *Vijla, *realGammaFai, *imagGammaFai;
  void sliceTensors();
  Tensor<double> &getSinglesContribution(const Map<3> &);
  Tensor<double> &getDoublesContribution(const Map<3> &);
  Tensor<double> &getEnergyDenominator(const Map<3> &);
};
} // namespace sisi4s

#endif
