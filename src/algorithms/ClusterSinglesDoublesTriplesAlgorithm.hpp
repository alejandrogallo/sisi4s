/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_TRIPLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_TRIPLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/FockVector.hpp>
#include <tcc/DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <ctf.hpp>

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
  class ClusterSinglesDoublesTriplesAlgorithm:
    public ClusterSinglesDoublesAlgorithm {

  public:

    ClusterSinglesDoublesTriplesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterSinglesDoublesTriplesAlgorithm();
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

    /**
     * \brief Computes and returns the residuum of the given amplitudes
     **/
    virtual PTR(FockVector<double>) getResiduum(
      const int iteration, const PTR(const FockVector<double>) &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes
     **/
    virtual PTR(FockVector<complex>) getResiduum(
      const int iteration, const PTR(const FockVector<complex>) &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the energy of the given amplitudes.
     **/
    //template <typename F>
    //F getEnergy(const PTR(const FockVector<F>) &amplitdues);

    /**
     * \brief Calculates an improved estimate of the amplitudes provided
     * the given the residuum.
     * \f$T_{ij\ldots}^{ab\ldots} = \frac{R_{ij\ldots}^{a\ldots}}
       {-\Delta_{ij\ldots}^{ab\ldots}}\f$
     * with \f$\Delta_{ij\ldots}^{ab\ldots} =
       \varepsilon_i+\ldots-\varepsilon_a-\ldots\f$.
     * \param[inout] residuum Fock vector, overwritten with new amplitudes.
     * \param[in] amplitudes Fock vector, previous amplitudes
     **/
    template <typename F>
    void estimateAmplitudesFromResiduum(
      const PTR(FockVector<F>) &residuum,
      const PTR(const FockVector<F>) &amplitudes
    );

    /**
     * \brief Calculates eps_a+eps_b+...-eps_i-eps_j-... into D^ab..._ij...
     **/
    template <typename F>
    void calculateExcitationEnergies(
      CTF::Tensor<F> &D, const std::string &indices
    );

    /**
     * \brief Dry run for amplitudesFromResiduum.
     * \param[in] R residuum tensor.
     **/
    template <typename F>
    void dryAmplitudesFromResiduum(cc4s::DryTensor<F> &R);

    template <typename F>
    PTR(FockVector<F>) createAmplitudes(
      std::initializer_list<std::string> amplitudeNames,
      std::initializer_list<std::initializer_list<int>> amplitudeLens,
      std::initializer_list<std::string> amplitudeIndices
    );

    template <typename F>
    void storeAmplitudes(
      const PTR(const FockVector<F>) &amplitudes,
      std::initializer_list<std::string> names
    );

    /**
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    std::string getDataName(const std::string &type, const std::string &data);
  };
}

#endif


