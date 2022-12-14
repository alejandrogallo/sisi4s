#ifndef PERTURBATIVE_TRIPLES_DEFINED
#define PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
class PerturbativeTriples : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(PerturbativeTriples);
  PerturbativeTriples(std::vector<Argument> const &argumentList);
  virtual ~PerturbativeTriples();
  /**
   * \brief Calculates perturbative triples correction. Routine based on
   * Helgaker book.
   */
  virtual void run();

  virtual void runInMemory();

  /**
   * \brief Calculates perturbative triples correction. Routine based on Piecuch
   * paper.
   */
  virtual void runPiecuch();
  /**
   * \brief Dry run for perturbative triples correction based on Helgaker book.
   */
  virtual void dryRun();
  /**
   * \brief Dry run for perturbative triples correction based on Piecuch paper.
   */
  virtual void dryRunPiecuch();
};
} // namespace sisi4s

#endif
