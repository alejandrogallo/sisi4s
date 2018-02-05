#ifndef EIGEN_SYSTEM_DAVIDSON_DEFINED
#define EIGEN_SYSTEM_DAVIDSON_DEFINED

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <math/MathFunctions.hpp>

#include <vector>
#include <iomanip>
#include <utility>
#include <algorithm>

namespace cc4s {
  template <typename H, typename P, typename V>
  class EigenSystemDavidson {
  public:
    typedef typename V::FieldType F;

    /**
     * \brief ...
     * \param[in] h object representing the matrix whose eigen system is sought
     * offering the following method:
     * V rightApply(V &v);
     * Also, if the dual version of the Davidson algorithm is to be used
     * h should offer the leftApply method
     * V leftApply(V &v);
     * returning the action of h on v standing right of h.
     * \param[in] eigenVectorsCount the number of eigenvalues and vectors to be
     * computed.
     * \param[in] p object representing the preconditioner
     * offering the following method:
     * std::vector<V> getInitialBasis(int eigenVectorsCount);
     * returning an initial guess the basis consisting of eigenVectorsCount
     * vectors.
     * V getCorrection(const complex eigenValue, V &residuum);
     * computing the estimated correction for the k-th eigenvector given its
     * eigenvalue and the residuum of the current k-th estimated eigenvector.
     * \param[in] tolerance the targeted relative tolerance to be met
     * by all residua.
     * \param[in] maxBasisSize the maximum allowed number of vectors
     * representing eigenvectors.
     * \param[in] dualVersion If the dual version of the algorithm is to be
     * used. The dual version of the algorithm calculates both the right and
     * left eigenvectors of the h object. If the dual version is not to be used
     * then only the right eigenvectors are to be calculated.
     * \param[in] refreshIterations This vector of integers represents
     * the iterations where the refreshment should be done. Refreshment meaning
     * that the trial Davidson basis should be thrown away and only keep
     * the current approximations to the eigenvectors.
     **/
    EigenSystemDavidson(
      H *h_,
      const int eigenVectorsCount_,
      P *p_,
      const double tolerance_,
      const unsigned int maxBasisSize_,
      const unsigned int maxIterations_,
      const unsigned int minIterations_
    ):
      h(h_),
      eigenVectorsCount(eigenVectorsCount_),
      p(p_),
      tolerance(tolerance_),
      maxBasisSize(maxBasisSize_),
      maxIterations(maxIterations_),
      minIterations(minIterations_),
      eigenValues(eigenVectorsCount_)
    {
    }

    virtual void run() = 0;

    const std::vector<complex> &getEigenValues() const {
      return eigenValues;
    }

    const std::vector<V> &getRightEigenVectors() const {
      return rightEigenVectors;
    }

    const std::vector<V> &getLeftEigenVectors() const {
      return leftEigenVectors;
    }

    /**
     * \brief Controls if the Davidson algorithm should make a refreshment
     * of the basis whenever the maximal basis size has been attained.
     * \param[in] value If value is true, it will be refreshed, otherwise not.
     */
    void refreshOnMaxBasisSize(const bool value) {
      refreshOnMaxBasisSizeValue = value;
    }

    /**
     * \brief Check wether or not the basis should be refreshed whenever
     * the maximal basis size has been reached.
     */
    bool refreshOnMaxBasisSize() {
      return refreshOnMaxBasisSizeValue;
    }


  protected:
    H *h;
    int eigenVectorsCount;
    P *p;
    double tolerance = 1E-14;
    unsigned int maxBasisSize = 1000;
    unsigned int maxIterations = 1000;
    unsigned int minIterations = 1;
    std::vector<int> refreshIterations = std::vector<int>{{}};
    bool refreshOnMaxBasisSizeValue = false;
    std::vector<complex> eigenValues;
    std::vector<V> rightEigenVectors;
    std::vector<V> leftEigenVectors;

  };



template <typename H, typename P, typename V>
class EigenSystemDavidsonMono: public EigenSystemDavidson<H,P,V> {
  public:
    typedef typename V::FieldType F;


    EigenSystemDavidsonMono(
      H *h_,
      const int eigenVectorsCount_,
      P *p_,
      const double tolerance_,
      const unsigned int maxBasisSize_,
      const unsigned int maxIterations_,
      const unsigned int minIterations_
    ): EigenSystemDavidson<H,P,V>(
      h_,
      eigenVectorsCount_,
      p_,
      tolerance_,
      maxBasisSize_,
      maxIterations_,
      minIterations_
    ) {}

    void run() {
      LOG(1,"Davidson").flags(
        std::ios::right | std::ios::scientific | std::ios::showpos
      );
      // get inital estimates for rEV = initial B matrix
      std::cout << this->maxBasisSize << std::endl;
      this->rightEigenVectors = this->p->getInitialBasis(this->eigenVectorsCount);

      std::vector<V> rightBasis( this->rightEigenVectors );

      // begin convergence loop
      double rms;
      unsigned int iterationCount(0);
      do {
        LOG(1,"Davidson") << "iteration=" << (iterationCount+1) << std::endl;

        // Check if a refreshment should be done
        if (
          std::find(
            this->refreshIterations.begin(),
            this->refreshIterations.end(),
            iterationCount + 1
          ) != this->refreshIterations.end()
          ||
          ( this->refreshOnMaxBasisSize() &&
            rightBasis.size() >= this->maxBasisSize )
        ) {
          LOG(1,"Davidson") << "Refreshing BASIS" << std::endl;
          rightBasis.resize(this->eigenVectorsCount);
          this->eigenValues.resize(this->eigenVectorsCount);
          for (unsigned int i(0) ; i < rightBasis.size() ; i++) {
            rightBasis[i] = this->rightEigenVectors[i];
          }
        }

        // compute reduced H by projection onto subspace spanned by rightBasis
        LapackMatrix<complex> reducedH(rightBasis.size(), rightBasis.size());

        for (unsigned int j(0); j < rightBasis.size(); ++j) {
          V HBj( this->h->rightApply(rightBasis[j]) );
          for (unsigned int i(0); i < rightBasis.size(); ++i) {
            //V HBi( h->rightApply(rightBasis[i]) );
            //reducedH(i,j) = HBi.dot(HBj);
            reducedH(i,j) = rightBasis[i].dot(HBj);
          }
        }

        // compute K lowest reduced eigenvalues and vectors of reduced H
        LapackMatrix<complex> reducedEigenVectors(
          rightBasis.size(), rightBasis.size()
        );
        LapackGeneralEigenSystem<complex> reducedEigenSystem(reducedH);

        // begin rightBasis extension loop for each k
        rms = 0.0;
        for (unsigned int k(0); k < this->eigenValues.size(); ++k) {
          // get estimated eigenvalue
          this->eigenValues[k] = reducedEigenSystem.getEigenValues()[k];

          // compute estimated eigenvector by expansion in rightBasis
          this->rightEigenVectors[k] *= F(0);
          for (int b(0); b < reducedH.getColumns(); ++b) {
            this->rightEigenVectors[k] +=
              rightBasis[b] * ComplexTraits<F>::convert(
                reducedEigenSystem.getRightEigenVectors()(b,k)
              );
          }
          double rightNorm(
            this->rightEigenVectors[k].dot(this->rightEigenVectors[k])
          );

          LOG(1,"Davidson") << "Right norm [" << k << "] "
                            << rightNorm << std::endl;
          LOG(1,"Davidson") << "EV         [" << k << "] "
                            << std::setprecision(15) << std::setw(23)
                            << this->eigenValues[k] << std::endl;



          // compute residuum
          V residuum( this->h->rightApply(this->rightEigenVectors[k]) );
          residuum -=
            this->rightEigenVectors[k] * ComplexTraits<F>::convert(
              //std::sqrt(eigenValues[k])
              this->eigenValues[k]
            );
          rms += std::real(residuum.dot(residuum)) /
            std::real(this->rightEigenVectors[k].dot(this->rightEigenVectors[k]));

          // compute correction using preconditioner
          V correction( this->p->getCorrection(this->eigenValues[k], residuum) );

          // orthonormalize and append to rightBasis
          for (unsigned int b(0); b < rightBasis.size(); ++b) {
            correction -= rightBasis[b] * rightBasis[b].dot(correction);
          }
          typename V::FieldType correction_norm(
            std::sqrt(correction.dot(correction))
          );
          if (std::abs(correction_norm) < 1E-6) continue;
          correction *= 1 / correction_norm;
          rightBasis.push_back(correction);
          LOG(1,"Davidson") << "Basis size " << rightBasis.size() << std::endl;
        }
        ++iterationCount;
        // end rightBasis extension loop
      } while (
        iterationCount+1 <= this->minIterations    || (
          rms >= this->eigenVectorsCount * this->tolerance &&
          rightBasis.size() <= this->maxBasisSize    &&
          iterationCount+1 <= this->maxIterations
        )
      );
      // end convergence loop
      if (rightBasis.size() > this->maxBasisSize) {
        //throw EXCEPTION("Failed to reach convergence");
      }
    }


};



} // namespace cc4s

#endif

/*
  MySimpleVector v;
  MySimpleMatrix m;
  MySimplePreconditioner<MySimpleMatrix> p
  EigenSystemDavidson<MyVectorType> eigenSystem(m, 4, p);
  eigenSystem.getEigenValues()[0];
  eigenSystem.getRightEigenVectors()[0];
*/
