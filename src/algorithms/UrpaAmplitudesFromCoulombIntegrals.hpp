#ifndef URPA_ALGORITHM__DEFINED
#define URPA_ALGORITHM__DEFINED
#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

namespace sisi4s {
class UrpaAmplitudesFromCoulombIntegrals
    : public ClusterSinglesDoublesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UrpaAmplitudesFromCoulombIntegrals);

  UrpaAmplitudesFromCoulombIntegrals(std::vector<Argument> const &argumentList)
      : ClusterSinglesDoublesAlgorithm(argumentList) {}
  ~UrpaAmplitudesFromCoulombIntegrals(){};

  std::string getAbbreviation() override { return "Urpa"; }

protected:
  /**
   * \brief Implements the iterate method with the DRCCD iteration.
   * \param[in] i Iteration number
   */
  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes);

  virtual PTR(FockVector<complex>)
  getResiduum(const int iteration,
              const PTR(const FockVector<complex>) &amplitudes);

  template <typename F>
  PTR(FockVector<F>)
  getResiduumTemplate(const int iteration,
                      const PTR(const FockVector<F>) &amplitudes);
};

} // namespace sisi4s

#endif
