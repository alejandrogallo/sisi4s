#ifndef UCCSDT_FROM_COULOMB_INTEGRALS_DEFINED
#define UCCSDT_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

namespace sisi4s {
class UccsdtAmplitudesFromCoulombIntegrals
    : public ClusterSinglesDoublesTriplesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UccsdtAmplitudesFromCoulombIntegrals);
  UccsdtAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList);
  virtual ~UccsdtAmplitudesFromCoulombIntegrals();

  virtual void run();
  virtual std::string getAbbreviation() { return "Uccsdt"; }

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

  template <typename F>
  PTR(FockVector<F>) getResiduumSth(const int iteration,
                                    const PTR(const FockVector<F>) &amplitudes);
};

} // namespace sisi4s

#endif
