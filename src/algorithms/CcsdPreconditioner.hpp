#ifndef CCSD_PRE_CONDITIONER_DEFINED
#define CCSD_PRE_CONDITIONER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>
#include <util/SharedPointer.hpp>

namespace cc4s {

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   */
  template <typename F>
    class CcsdPreconditioner {
    public:
      typedef SDFockVector<F> V;

      CcsdPreconditioner(){}
      ~CcsdPreconditioner(){}

      /**
       * \brief Setters for the main tensors
       */
      CcsdPreconditioner&
      setTai(CTF::Tensor<F> *t) { Tai = t; return *this;}
      CcsdPreconditioner&
      setTabij(CTF::Tensor<F> *t) { Tabij = t; return *this;}
      CcsdPreconditioner&
      setFij(CTF::Tensor<F> *t) { Fij = t; return *this;}
      CcsdPreconditioner&
      setFab(CTF::Tensor<F> *t) { Fab = t; return *this;}
      CcsdPreconditioner&
      setVabcd(CTF::Tensor<F> *t) { Vabcd = t; return *this;}
      CcsdPreconditioner&
      setViajb(CTF::Tensor<F> *t) { Viajb = t; return *this;}
      CcsdPreconditioner&
      setVijab(CTF::Tensor<F> *t) { Vijab = t; return *this;}
      CcsdPreconditioner&
      setVijkl(CTF::Tensor<F> *t) { Vijkl = t; return *this;}

      CcsdPreconditioner&
      setRandom(bool t) { preconditionerRandom = t; return *this;}
      CcsdPreconditioner&
      setRandomSigma(double t) { preconditionerRandomSigma=t; return *this;}

      /**
       * \brief Get initial basis
       * \param[in] eigenVectorsCount Number of eigen vectors
       */
      std::vector<V> getInitialBasis(int eigenVectorsCount);

      V getCorrection(const complex eigenValue, V &residuum);

      SDTFockVector<F>
      getCorrection(const complex eigenValue, SDTFockVector<F> &residuum);

      void calculateDiagonal();
      PTR(V) getDiagonal() {
        if (!diagonalH)
         calculateDiagonal();
        return diagonalH;
      }

    private:
      PTR(V) diagonalH;
      CTF::Tensor<F> *Fij;
      CTF::Tensor<F> *Fab;
      CTF::Tensor<F> *Tai = nullptr;
      CTF::Tensor<F> *Tabij = nullptr;
      CTF::Tensor<F> *Vabcd = nullptr;
      CTF::Tensor<F> *Viajb = nullptr;
      CTF::Tensor<F> *Vijab = nullptr;
      CTF::Tensor<F> *Vijkl = nullptr;

      /**
       * Wether or not to use random preconditioners.
       */
      bool preconditionerRandom = false;

      /**
       * The standard deviation used in the normal distribution to create
       * random preconditioners.
       */
      double preconditionerRandomSigma = 1.0;

  };

}

#endif


