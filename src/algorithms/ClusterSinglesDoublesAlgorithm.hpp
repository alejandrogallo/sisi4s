#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <util/Tensor.hpp>

#include <string>
#include <vector>

namespace sisi4s {

#define CLUSTER_SINGLES_DOUBLES_INSPEC                                         \
  {"amplitudesConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-5)},        \
      {"energyConvergence", SPEC_VALUE_DEF("TODO: DOC", double, 1e-6)},        \
      {"levelShift", SPEC_VALUE_DEF("TODO: DOC", double, 0.0)},                \
      {"antisymmetrize", SPEC_VALUE_DEF("TODO: DOC", bool, false)},            \
      {"integralsSliceSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},        \
      {"maxIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 16)},             \
      {"unrestricted", SPEC_VALUE_DEF("TODO: DOC", bool, false)},              \
      {"mixer", SPEC_VALUE_DEF("TODO: DOC", std::string, "LinearMixer")},      \
      {"maxResidua", SPEC_VALUE_DEF("TODO: DOC", int64_t, 4)},                 \
      {"mixingRatio", SPEC_VALUE_DEF("TODO: DOC", double, 1.0)},               \
      {"distinguishable", SPEC_VALUE_DEF("TODO: DOC", bool, false)},           \
      {"OnlyPPL", SPEC_VALUE_DEF("TODO: DOC", bool, false)},                   \
      {"PPL", SPEC_VALUE_DEF("TODO: DOC", bool, true)},                        \
      {"CoulombFactors", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},  \
      {"CoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},   \
      {"FactorOrbitals", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},  \
      {"OutgoingFactorOrbitals",                                               \
       SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},                    \
      {"HoleEigenEnergies",                                                    \
       SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},                  \
      {"ParticleEigenEnergies",                                                \
       SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},                  \
      {"HHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<F> *)},          \
      {"HPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<F> *)},                  \
      {"initialSinglesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<F> *)},      \
      {"initialDoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<F> *)},      \
      {"PPHHCoulombIntegrals",                                                 \
       SPEC_VARIN("TODO: DOC", Tensor<F> *)->require()},                       \
  {                                                                            \
    "CoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)        \
  }

#define CLUSTER_SINGLES_DOUBLES_OUTSPEC                                        \
  {"SinglesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)}, {              \
    "DoublesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)                 \
  }

/**
 * \brief Contains all the necessary tools for an algorithm with
 * singles and doubles amplitudes. It calculates the energy from the amplitudes
 * \f$T_{a}^{i}\f$ and \f$T_{ab}^{ij}\f$ and the Coulomb integrals
 *\f$V_{ij}^{ab}\f$. For calculating the amplitudes it calls the iteration
 *routine of the actual algorithm.
 **/
class ClusterSinglesDoublesAlgorithm : public Algorithm {
public:
  using Algorithm::Algorithm;

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
  static double constexpr DEFAULT_ENERGY_CONVERGENCE = 1E-6;
  static double constexpr DEFAULT_AMPLITUDES_CONVERGENCE = 1E-5;

  static double constexpr DEFAULT_LEVEL_SHIFT = 0.0;

protected:
  template <typename F>
  F run();

  /**
   * \brief Computes and returns the residuum of the given amplitudes
   **/
  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes) = 0;

  /**
   * \brief Computes and returns the residuum of the given amplitudes
   **/
  virtual PTR(FockVector<complex>)
  getResiduum(const int iteration,
              const PTR(const FockVector<complex>) &amplitudes) = 0;

  /**
   * \brief Computes and returns the energy of the given amplitudes.
   **/
  template <typename F>
  F getEnergy(const PTR(const FockVector<F>) &amplitdues);

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
  void
  estimateAmplitudesFromResiduum(const PTR(FockVector<F>) &residuum,
                                 const PTR(const FockVector<F>) &amplitudes);

  /**
   * \brief Calculates eps_a+eps_b+...-eps_i-eps_j-... into D^ab..._ij...
   **/
  template <typename F>
  void calculateExcitationEnergies(Tensor<F> &D, const std::string &indices);

  /**
   * \brief Dry run for amplitudesFromResiduum.
   * \param[in] R residuum tensor.
   **/
  template <typename F>
  void dryAmplitudesFromResiduum(sisi4s::DryTensor<F> &R);

  template <typename F>
  PTR(FockVector<F>)
  createAmplitudes(std::vector<std::string> amplitudeNames,
                   std::vector<std::vector<TensorIndex>> amplitudeLens,
                   std::vector<std::string> amplitudeIndices);

  template <typename F>
  void storeAmplitudes(const PTR(const FockVector<F>) &amplitudes,
                       std::vector<std::string> names);

  /**
   * \brief Calculates and returns one slice Xxycd of the Coulomb integrals
   * \f$V_{cd}^{ab}\f$ coupled to the singles amplitudes. The indices x and y
   * are restricted to the range {No+a, ..., No+a+No-1} and {No+b, ...,
   * No+b+No-1}, respectively. The caller is responsible for deleting the
   * dynamically allocated result tensor. \param[in] a 1st sliced dimension (x).
   * \param[in] b 2nd sliced dimension (y).
   * \param[in] integralsSliceSize slicing rank.
   * \param[out] Xxycd sliced coupled Coulomb integrals Xabcd
   */
  Tensor<double> *
  sliceCoupledCoulombIntegrals(const PTR(const FockVector<double>) &amplitudes,
                               int a,
                               int b,
                               int integralsSliceSize);

  Tensor<sisi4s::complex> *
  sliceCoupledCoulombIntegrals(const PTR(const FockVector<complex>) &amplitudes,
                               int a,
                               int b,
                               int integralsSliceSize);

  /**
   * \brief Calculates and returns one slice Fabij of the residuum
   * from the dressed Coulomb factors. The slice is computed from
   * Rx and Ry and are restricted to the
   * range {a, ..., factorsSliceSize+a-1} and {b, ..., factorsSliceSize+b-1},
   * respectively. The caller is responsible for deleting the dynamically
   * allocated result tensor. \param[in] a 1st sliced dimension (Rx). \param[in]
   * b 2nd sliced dimension (Ry). \param[in] factorsSliceSize slicing rank of
   * NR. \param[out] Fabij sliced Residuum
   */
  Tensor<double> *sliceAmplitudesFromCoupledCoulombFactors(
      const PTR(const FockVector<double>) &amplitudes,
      int a,
      int b,
      int factorsSliceSize);
  Tensor<sisi4s::complex> *sliceAmplitudesFromCoupledCoulombFactors(
      const PTR(const FockVector<complex>) &amplitudes,
      int a,
      int b,
      int factorsSliceSize);

  /**
   * \brief Adds the given slice of the residuum tensor Rxyij to the
   * entire residuum tensor Rabij at the respective index range.
   * \param[in] a0 1st sliced dimension (x).
   * \param[in] b0 2nd sliced dimension (y).
   * \param[in] Rxyij sliced residuum
   * \param[in] Rabij entire residuum.
   */
  template <typename F>
  void sliceIntoResiduum(Tensor<F> &Rxyij, int a0, int b0, Tensor<F> &Rabij);

  /**
   * \brief The abbreviation of the algorithm in capital letters.
   **/
  std::string getCapitalizedAbbreviation();

  std::string getDataName(const std::string &type, const std::string &data);
};
} // namespace sisi4s

#endif
