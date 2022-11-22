/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef UPERTURBATIVE_TRIPLES_DEFINED
#define UPERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
class UPerturbativeTriples : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UPerturbativeTriples);
  UPerturbativeTriples(std::vector<Argument> const &argumentList);
  virtual ~UPerturbativeTriples();
  /**
   * \brief Calculates perturbative triples correction. Routine based on
   * Helgaker book.
   */
  virtual void run();

  /**
   * \brief Dry run for perturbative triples correction based on Helgaker book.
   */
  virtual void dryRun();
};
} // namespace sisi4s

#endif
