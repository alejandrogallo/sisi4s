/*Copyright (c) 2022, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSDTENERGYFROMCOULOMBINTEGRALS_HPP_
#define CCSDTENERGYFROMCOULOMBINTEGRALS_HPP_

#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>

namespace sisi4s {
class CcsdtEnergyFromCoulombIntegrals
    : public ClusterSinglesDoublesTriplesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdtEnergyFromCoulombIntegrals);
  CcsdtEnergyFromCoulombIntegrals(std::vector<Argument> const &argumentList)
      : ClusterSinglesDoublesTriplesAlgorithm(argumentList){};
  virtual ~CcsdtEnergyFromCoulombIntegrals(){};

  virtual std::string getAbbreviation() { return "Ccsdt"; }

  static int64_t constexpr DEFAULT_SLICE_SIZE = -1;
  static int64_t constexpr DEFAULT_DISTINGUISHABLE = 0;

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
