#ifndef ITERATIVE_PSEUDO_INVERSE_DEFINED
#define ITERATIVE_PSEUDO_INVERSE_DEFINED

#include <math/Complex.hpp>
#include <DryTensor.hpp>
#include <util/Tensor.hpp>
#include <random>

namespace sisi4s {
  template <typename F>
  class IterativePseudoInverse {
  public:
    IterativePseudoInverse(Tensor<F> const &matrix, F accuracy=1e-10);
    Tensor<F> &get();

    static void test(CTF::World *world);
  protected:
    void iterate(F accuracy = 1e-10);
    void iterateQuadratically(F accuracy = 1e-10);
    static void setRandom(
      F &value,
      std::mt19937 &random, std::normal_distribution<F> &normalDistribution
    );
    static void generateHilbertMatrix(Tensor<F> &matrix);

    Tensor<F> matrix;
    Tensor<F> square, inverse;
    F alpha;
  };

  template <typename F>
  class DryIterativePseudoInverse {
  public:
    DryIterativePseudoInverse(DryTensor<F> const &matrix);
    DryTensor<F> &get();

  protected:
    DryTensor<F> matrix;
    DryTensor<F> square, inverse;
  };
}

#endif

