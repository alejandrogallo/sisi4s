/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RALS_FTOD_RANK_DECOMPOSITION_DEFINED
#define RALS_FTOD_RANK_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <math/Complex.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief This algorithm provides a tensor rank decomposition of the
   * Fourier tranformed overlap densities.
   */
  class RalsFtodRankDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RalsFtodRankDecomposition);
    RalsFtodRankDecomposition(std::vector<Argument> const &argumentList);
    virtual ~RalsFtodRankDecomposition();
    virtual void run();
      
    /**
     * \brief the rank of the tensor rank decomposition
     */
    int64_t rank;
    double R;
    CTF::Tensor<complex> *chi, *chi0;
    CTF::Matrix<complex> *x, *gamma;

    static void test(CTF::World *world);
  protected:
    void fit(double lambda);

    void fitAls(
      char const *indicesChi,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA
    );
    double fitRals(
      char const *indicesChi,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA,
      double lambda
    );

    void normalizeX();
    void realizeX();
  };
}

#endif

