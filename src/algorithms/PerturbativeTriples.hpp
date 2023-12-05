#ifndef PERTURBATIVE_TRIPLES_DEFINED
#define PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
DEFINE_ALGORITHM_HEADER(

    PerturbativeTriples,

    virtual void runInMemory();

    /**
     * \brief Calculates perturbative triples correction. Routine based on
     * Piecuch paper.
     */
    virtual void runPiecuch();
    /**
     * \brief Dry run for perturbative triples correction based on Piecuch
     * paper.
     */
    virtual void dryRunPiecuch(););
} // namespace sisi4s

#endif
