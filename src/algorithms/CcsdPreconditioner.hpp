#ifndef CCSD_PRE_CONDITIONER_DEFINED
#define CCSD_PRE_CONDITIONER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>
#include <math/Complex.hpp>
#include <util/SharedPointer.hpp>

namespace sisi4s {

/**
 * \brief Implements the diagonal preconditionar for the davidson method
 * \tparam F It is the field variable to be used, in general it will be
 * complex
 */
template <typename F>
class CcsdPreconditioner {
public:
  typedef SDFockVector<F> V;

  CcsdPreconditioner() {}
  ~CcsdPreconditioner() {}

  /**
   * \brief Setters for the main tensors
   */
  CcsdPreconditioner &setTai(Tensor<F> *t) {
    Tai = t;
    return *this;
  }
  CcsdPreconditioner &setTabij(Tensor<F> *t) {
    Tabij = t;
    return *this;
  }
  CcsdPreconditioner &setFij(Tensor<F> *t) {
    Fij = t;
    return *this;
  }
  CcsdPreconditioner &setFab(Tensor<F> *t) {
    Fab = t;
    return *this;
  }
  CcsdPreconditioner &setVabcd(Tensor<F> *t) {
    Vabcd = t;
    return *this;
  }
  CcsdPreconditioner &setViajb(Tensor<F> *t) {
    Viajb = t;
    return *this;
  }
  CcsdPreconditioner &setVijab(Tensor<F> *t) {
    Vijab = t;
    return *this;
  }
  CcsdPreconditioner &setVijkl(Tensor<F> *t) {
    Vijkl = t;
    return *this;
  }

  CcsdPreconditioner &setSpinFlip(bool t) {
    spinFlip = t;
    return *this;
  }

  CcsdPreconditioner &setRandom(bool t) {
    preconditionerRandom = t;
    return *this;
  }
  CcsdPreconditioner &setRandomSigma(double t) {
    preconditionerRandomSigma = t;
    return *this;
  }

  /**
   * \brief Get initial basis
   * \param[in] eigenVectorsCount Number of eigen vectors
   */
  std::vector<SDFockVector<F>> getInitialBasis(int eigenVectorsCount);

  SFockVector<F> getCorrection(const complex eigenValue,
                               SFockVector<F> &residuum);

  SDFockVector<F> getCorrection(const complex eigenValue,
                                SDFockVector<F> &residuum);

  SDTFockVector<F> getCorrection(const complex eigenValue,
                                 SDTFockVector<F> &residuum);

  void calculateDiagonal();
  PTR(V) getDiagonal() {
    if (!diagonalH) calculateDiagonal();
    return diagonalH;
  }

  PTR(SDFockVector<F>) diagonalH;
  Tensor<F> *Fij;
  Tensor<F> *Fab;
  Tensor<F> *Tai = nullptr;
  Tensor<F> *Tabij = nullptr;
  Tensor<F> *Vabcd = nullptr;
  Tensor<F> *Viajb = nullptr;
  Tensor<F> *Vijab = nullptr;
  Tensor<F> *Vijkl = nullptr;

private:
  /**
   * Wether or not to use random preconditioners.
   */
  bool preconditionerRandom = false;

  // wether or not to use spin flip in the filtering
  bool spinFlip = false;

  /**
   * The standard deviation used in the normal distribution to create
   * random preconditioners.
   */
  double preconditionerRandomSigma = 1.0;
};

template <typename F>
class IPCcsdPreconditioner : public CcsdPreconditioner<F> {
public:
  void calculateDiagonal();

  std::vector<SDFockVector<F>> getInitialBasis(int eigenVectorsCount);

  SDFockVector<F> getCorrection(const complex eigenValue,
                                SDFockVector<F> &residuum);
};

template <typename F>
class EACcsdPreconditioner : public CcsdPreconditioner<F> {
public:
  void calculateDiagonal();

  std::vector<SDFockVector<F>> getInitialBasis(int eigenVectorsCount);

  SDFockVector<F> getCorrection(const complex eigenValue,
                                SDFockVector<F> &residuum);
};

} // namespace sisi4s

#endif
