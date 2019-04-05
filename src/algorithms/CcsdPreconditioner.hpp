#ifndef CCSD_PRE_CONDITIONER_DEFINED
#define CCSD_PRE_CONDITIONER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>

namespace cc4s {

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   */
  template <typename F>
    class CcsdPreconditioner {
      public:
        typedef CcsdFockVector<F> V;

        /**
         * Wether or not to use random preconditioners.
         */
        bool preconditionerRandom = false;

        /**
         * The standard deviation used in the normal distribution to create
         * random preconditioners.
         */
        double preconditionerRandomSigma = 1.0;

        CcsdPreconditioner(){}

        /**
         * \brief Constructor for the preconditioner.
         */
        CcsdPreconditioner (
          CTF::Tensor<F> &Tai,
          CTF::Tensor<F> &Tabij,
          CTF::Tensor<F> &Fij,
          CTF::Tensor<F> &Fab,
          CTF::Tensor<F> &Vabcd,
          CTF::Tensor<F> &Viajb,
          CTF::Tensor<F> &Vijab,
          CTF::Tensor<F> &Vijkl
        );

        /**
         * \brief Get initial basis
         * \param[in] eigenVectorsCount Number of eigen vectors
         */
        std::vector<V> getInitialBasis(int eigenVectorsCount);

        V getCorrection(const complex eigenValue, V &residuum);

        CcsdtFockVector<F>
        getCorrection(const complex eigenValue, CcsdtFockVector<F> &residuum);

        V getDiagonalH() const { return diagonalH; }

  protected:
    V diagonalH;

  };


}

#endif


