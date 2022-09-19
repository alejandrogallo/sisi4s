/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REGULARIZED_ALTERNATING_LEAST_SQUARES_DEFINED
#define REGULARIZED_ALTERNATING_LEAST_SQUARES_DEFINED

#include <math/Complex.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
  class AlternatingLeastSquaresRegularizationEstimator {
  public:
    AlternatingLeastSquaresRegularizationEstimator(): swampingThreshold(0.0) {
    }
    AlternatingLeastSquaresRegularizationEstimator(
      double swampingThreshold_, double regularizationFriction_,
      double initialLambda_
    ):
      swampingThreshold(swampingThreshold_),
      regularizationFriction(regularizationFriction_),
      lambda(initialLambda_)
    { }
    virtual ~AlternatingLeastSquaresRegularizationEstimator() {
    }
    double getSwampingThreshold() {
      return swampingThreshold;
    }
    virtual double getLambda() {
      return lambda;
    }
    virtual void update(double const swampingFactor) {
      double s(swampingFactor / swampingThreshold);
      double estimatedLambda(lambda * s*s);
      lambda =
        (1-regularizationFriction)*estimatedLambda +
        regularizationFriction*lambda;
    }
  protected:
    double swampingThreshold, regularizationFriction;
    double lambda;
  };

  class NoRegularizationEstimator:
    public AlternatingLeastSquaresRegularizationEstimator
  {
  public:
    virtual ~NoRegularizationEstimator() {
    }
    virtual double getLambda() {
      return 0.0;
    }
    virtual void update(double const swampingFactor) {
    }
  };

  template <typename F=double>
  void fitAlternatingLeastSquaresFactor(
    Tensor<F> &T, char const *indicesT,
    Tensor<F> &B, char const idxB,
    Tensor<F> &C, char const idxC,
    Tensor<F> &A, char const idxA
  );

  template <typename F=double>
  void fitRegularizedAlternatingLeastSquaresFactor(
    Tensor<F> &T, char const *indicesT,
    Tensor<F> &B, char const idxB,
    Tensor<F> &C, char const idxC,
    Tensor<F> &A, char const idxA,
    AlternatingLeastSquaresRegularizationEstimator *regularizationEstimatorA
  );

  template <typename F=double>
  void dryFitRegularizedAlternatingLeastSquaresFactor(
    DryTensor<F> &T, char const *indicesT,
    DryTensor<F> &B, char const idxB,
    DryTensor<F> &C, char const idxC,
    DryTensor<F> &A, char const idxA
  );
}

#endif

