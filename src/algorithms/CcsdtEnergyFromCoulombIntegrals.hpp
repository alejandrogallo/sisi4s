#ifndef CCSDTENERGYFROMCOULOMBINTEGRALS_HPP_
#define CCSDTENERGYFROMCOULOMBINTEGRALS_HPP_

#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>

namespace sisi4s {
class CcsdtEnergyFromCoulombIntegrals
    : public ClusterSinglesDoublesTriplesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdtEnergyFromCoulombIntegrals);
  using ClusterSinglesDoublesTriplesAlgorithm::
      ClusterSinglesDoublesTriplesAlgorithm;
  static AlgorithmInputSpec spec;
  virtual std::string getAbbreviation() { return "Ccsdt"; }
  template <typename F>
  PTR(FockVector<F>)
  getResiduumTemplate(const int iterationStep,
                      const PTR(const FockVector<F>) &amplitudes);

  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes);
  virtual PTR(FockVector<complex>)
  getResiduum(const int iteration,
              const PTR(const FockVector<complex>) &amplitudes);
};
} // namespace sisi4s

#endif
