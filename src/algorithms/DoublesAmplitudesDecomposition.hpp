#ifndef DOUBLES_AMPLITUDES_DECOMPOSITION_DEFINED
#define DOUBLES_AMPLITUDES_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <memory>

namespace sisi4s {
class DoublesAmplitudesDecomposition : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(DoublesAmplitudesDecomposition);
  DoublesAmplitudesDecomposition(std::vector<Argument> const &argumentList);
  virtual ~DoublesAmplitudesDecomposition();
  /**
   * \brief Calculates left singular vectors of the particle-hole Coulomb vertex
   * \f$\tilde\Gamma^a_{iG}\f$.
   */
  virtual void run();
  /**
   * \brief Dry run for calculating the left singular vectors of
   * \f$\Gamma^a_{iG}\f$.
   */
  virtual void dryRun();

  static double constexpr DEFAULT_REDUCTION = 2.0;
  static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;

protected:
  void diagonlizeAmplitudes();
  void sliceLargestEigenValues();

  int Nv, No, NvNo, NF, lower, upper;
  double *lambdas;
  complex *sqrtLambdas;
  int64_t lambdasCount, *lambdaIndices;

  std::shared_ptr<Tensor<double>> Taibj, UaiF;
  std::shared_ptr<Tensor<complex>> LFai, sqrtLambdaF;
  Tensor<double> *LambdaF;
};
} // namespace sisi4s

#endif
