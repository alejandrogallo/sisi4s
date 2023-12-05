#ifndef CCSD_PERTURBATIVE_TRIPLES_DEFINED
#define CCSD_PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Permutation.hpp>
#include <util/SlicedCtfTensor.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction
 */
DEFINE_ALGORITHM_HEADER(

    CcsdPerturbativeTriples,

    int No,
    Nv;
    Tensor<double> * SVabc, *DVabc;
    Tensor<double> * realGammaFab, *imagGammaFab;
    SlicedCtfTensor<> * Tai, *Tabij, *Tabil;
    SlicedCtfTensor<> * Vabij, *Vijla, *realGammaFai, *imagGammaFai;
    void sliceTensors();
    Tensor<double> & getSinglesContribution(const Map<3> &);
    Tensor<double> & getDoublesContribution(const Map<3> &);
    Tensor<double> & getEnergyDenominator(const Map<3> &););
} // namespace sisi4s

#endif
