/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/FockVector.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <util/CTF.hpp>

#include <string>
#include <initializer_list>

namespace cc4s {
  /**
   * \brief Contains all the necessary tools for an algorithm with
   * singles, doubles and triples amplitudes.
   * It calculates the energy from the amplitudes
   * \f$T_{a}^{i}\f$, \f$T_{ab}^{ij}\f$ \f$T_{abc}^{ijk}\f$ and the Coulomb
   * integrals \f$V_{ij}^{ab}\f$.
   * For calculating the amplitudes it calls the iteration
   * routine of the actual algorithm.
   **/
  class ClusterSinglesDoublesTriplesQuadruplesAlgorithm:
    public ClusterSinglesDoublesAlgorithm {

  public:

    ClusterSinglesDoublesTriplesQuadruplesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterSinglesDoublesTriplesQuadruplesAlgorithm();
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

    /**
     * \brief Defines the default number of iterations (16).
     */
    static int constexpr DEFAULT_MAX_ITERATIONS = 16;

    static double constexpr DEFAULT_LEVEL_SHIFT = 0.0;

  protected:
    template <typename F>
    F run();

  };
}

#endif


