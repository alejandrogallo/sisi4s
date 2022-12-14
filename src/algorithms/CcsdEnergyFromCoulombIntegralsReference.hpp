#ifndef CCSD_ENERGY_FROM_COULOMB_INTEGRALS_REFERENCE_DEFINED
#define CCSD_ENERGY_FROM_COULOMB_INTEGRALS_REFERENCE_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

namespace sisi4s {
// this algorithm is now based on the ClusterSinglesDoublesAlgorithm
// inheriting its iteration and slicing functionality.
// Only the abstract (left out) methods getAbbreviation and iterate have
// to be implemented.
/**
 * \brief Implements the iteration routine for the Ccsd method. Calculates the
 * amplitudes \f$T_{a}^{i}\f$ and \f$T_{ab}^{ij}\f$ from the Coulomb
 * integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
 * V_{kl}^{ij}, V_{ka}^{ij}, V_{ci}^{ab}\f$ and \f$V_{cd}^{ab}\f$ (if given,
 * else slicing and the Coulomb Vertex \f$\Gamma_{pG}^q\f$  is used).
 */
class CcsdEnergyFromCoulombIntegralsReference
    : public ClusterSinglesDoublesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(CcsdEnergyFromCoulombIntegralsReference);
  CcsdEnergyFromCoulombIntegralsReference(
      std::vector<Argument> const &argumentList);
  virtual ~CcsdEnergyFromCoulombIntegralsReference();

  /**
   * \brief Returns the abbreviation of the routine (CCSD).
   * \return abbreviation of the routine
   */
  virtual std::string getAbbreviation() { return "Ccsd"; }

  static int64_t constexpr DEFAULT_SLICE_SIZE = -1;
  static int64_t constexpr DEFAULT_DISTINGUISHABLE = 0;

protected:
  /**
   * \brief Implements the iterate method with the CCSD iteration. Iteration
   * routine taken from So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001).
   * \param[in] i Iteration number
   */
  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes);
  virtual PTR(FockVector<complex>)
  getResiduum(const int iteration,
              const PTR(const FockVector<complex>) &amplitudes);
};
} // namespace sisi4s

#endif
