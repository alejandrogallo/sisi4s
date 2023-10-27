#ifndef CCSD_PRE_CONDITIONER_DEFINED
#define CCSD_PRE_CONDITIONER_DEFINED

#define _DEFINE_SETTER(type, name, default)                                    \
  CcsdPreconditioner &set##name(type t) {                                      \
    name = t;                                                                  \
    return *this;                                                              \
  }                                                                            \
  type name = default

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
  _DEFINE_SETTER(Tensor<F> *, Tai, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Tabij, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Fij, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Fab, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Vabcd, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Vijab, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Viajb, nullptr);
  _DEFINE_SETTER(Tensor<F> *, Vijkl, nullptr);

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

template <typename F>
class CISPreconditioner : public CcsdPreconditioner<F> {
public:
  using V = SFockVector<F>;

  PTR(V) diagonalH;
  void calculateDiagonal();

  std::vector<V> getInitialBasis(int eigenVectorsCount);

  V getCorrection(const complex eigenValue, V &residuum);
};

} // namespace sisi4s

#undef _DEFINE_SETTER

#endif
