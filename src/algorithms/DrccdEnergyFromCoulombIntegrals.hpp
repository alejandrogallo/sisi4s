#ifndef DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED
#define DRCCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

#include <util/SharedPointer.hpp>

namespace sisi4s {
// this algorithm is now based on the ClusterSinglesDoublesAlgorithm
// inheriting its iteration and slicing functionality.
// Only the abstract (left out) methods getAbbreviation and iterate have
// to be implemented.
/**
 * \brief Implements the iteration routine for the Drccd method. Calculates the
 * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
 * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
 */
class DrccdEnergyFromCoulombIntegrals : public ClusterSinglesDoublesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(DrccdEnergyFromCoulombIntegrals);
  DrccdEnergyFromCoulombIntegrals(std::vector<Argument> const &argumentList);
  virtual ~DrccdEnergyFromCoulombIntegrals();
  /**
   * \brief Returns the abbreviation of the routine (DRCCD).
   * \return abbreviation of the routine
   */
  virtual std::string getAbbreviation() { return "Drccd"; }

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
  PTR(FockVector<F>) getResiduum(const int iteration,
                                 const PTR(const FockVector<F>) &amplitudes);
};
} // namespace sisi4s

#endif
