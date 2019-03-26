/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all
 * rights reserved.*/
#ifndef CCSD_EQUATION_OF_MOTION_DAVIDSON
#define CCSD_EQUATION_OF_MOTION_DAVIDSON

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
  template <typename F = complex>
    class CcsdPreConditioner {
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


        /**
         * \brief Constructor for the preconditioner.
         */
        CcsdPreConditioner (
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

        V getDiagonalH() const { return diagonalH; }

  protected:
    V diagonalH;

  };

  class CcsdEquationOfMotionDavidson: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdEquationOfMotionDavidson);
    CcsdEquationOfMotionDavidson(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdEquationOfMotionDavidson();

    virtual void run();

    template<typename F>
    void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };
}

#endif

