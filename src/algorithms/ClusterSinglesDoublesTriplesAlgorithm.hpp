#ifndef CLUSTER_SINGLES_DOUBLES_TRIPLES_ALGORITHM_DEFINED
#define CLUSTER_SINGLES_DOUBLES_TRIPLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/FockVector.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <util/Tensor.hpp>

#include <string>
#include <initializer_list>

namespace sisi4s {
/**
 * \brief Contains all the necessary tools for an algorithm with
 * singles, doubles and triples amplitudes.
 * It calculates the energy from the amplitudes
 * \f$T_{a}^{i}\f$, \f$T_{ab}^{ij}\f$ \f$T_{abc}^{ijk}\f$ and the Coulomb
 * integrals \f$V_{ij}^{ab}\f$.
 * For calculating the amplitudes it calls the iteration
 * routine of the actual algorithm.
 **/
class ClusterSinglesDoublesTriplesAlgorithm
    : public ClusterSinglesDoublesAlgorithm {

public:
  using ClusterSinglesDoublesAlgorithm::ClusterSinglesDoublesAlgorithm;
  virtual ~ClusterSinglesDoublesTriplesAlgorithm() {}
  /**
   * \brief Calculates the energy of a ClusterSinglesDoubles algorithm
   */
  virtual void run();

  // TODO: dryRun

  /**
   * \brief Returns the abbreviation of the concrete algorithm, e.g. "Ccd",
   * "Dcd".
   */
  virtual std::string getAbbreviation() = 0;

  static double constexpr DEFAULT_LEVEL_SHIFT = 0.0;

protected:
  template <typename F>
  F run();
};
} // namespace sisi4s

#endif
