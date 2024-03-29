#include <math/ComplexTensor.hpp>

#include <util/Exception.hpp>
#include <util/Tensor.hpp>
#include <complex>

#define AssertCompatibleTensorShapes(c, r, i)                                  \
  {                                                                            \
    Assert((c).order == (r).order && (c).order == (i).order,                   \
           "Incompatible tensor orders");                                      \
    for (int k(0); k < c.order; ++k)                                           \
      Assert((c).lens[k] == (r).lens[k] && (c).lens[k] == (i).lens[k],         \
             "Incompatible tensor shapes");                                    \
  }

#define AssertCompatibleTensorShape(c, r)                                      \
  {                                                                            \
    Assert((c).order == (r).order, "Incompatible tensor orders");              \
    for (int k(0); k < c.order; ++k)                                           \
      Assert((c).lens[k] == (r).lens[k], "Incompatible tensor shapes");        \
  }

#define Indices(Tensor)                                                        \
  char indices[(Tensor).order];                                                \
  for (int i(0); i < (Tensor).order; ++i) { indices[i] = 'a' + i; }            \
  indices[(Tensor).order] = 0

void sisi4s::fromComplexTensor(Tensor<complex> &C,
                               Tensor<double> &R,
                               Tensor<double> &I) {
  AssertCompatibleTensorShapes(C, R, I);
  Indices(C);
  fromComplexTensor(C, R);
  I[indices] = CTF::Function<complex, double>(std::function<double(complex)>(
      [](complex c) { return c.imag(); }))(C[indices]);
}

void sisi4s::fromComplexTensor(Tensor<double> &C,
                               Tensor<double> &R,
                               Tensor<double> &I) {
  AssertCompatibleTensorShapes(C, R, I);
  Indices(C);
  R[indices] = C[indices];
  I[indices] = 0.0;
}

void sisi4s::fromComplexTensor(Tensor<complex> &C, Tensor<double> &R) {
  AssertCompatibleTensorShape(C, R);
  Indices(C);
  R[indices] = CTF::Function<complex, double>(std::function<double(complex)>(
      [](complex c) { return c.real(); }))(C[indices]);
}

void sisi4s::toComplexTensor(Tensor<double> &R,
                             Tensor<double> &I,
                             Tensor<complex> &C) {
  AssertCompatibleTensorShapes(C, R, I);
  Indices(C);
  toComplexTensor(R, C);
  CTF::Transform<double, complex>(
      std::function<void(double, complex &)>([](double i, complex &c) {
        // #ifdef INTEL_COMPILER
        //         c.imag() = i;
        // #else
        c.imag(i);
        // #endif
      }))(I[indices], C[indices]);
}

void sisi4s::toComplexTensor(Tensor<double> &R, Tensor<complex> &C) {
  AssertCompatibleTensorShape(C, R);
  Indices(C);
  CTF::Transform<double, complex>(
      std::function<void(double, complex &)>([](double r, complex &c) {
        // #ifdef INTEL_COMPILER
        //         c.real() = r;
        // #else
        c.real(r);
        c.imag(0);
        // #endif
      }))(R[indices], C[indices]);
}

void sisi4s::toComplexTensor(Tensor<double> &R, Tensor<double> &C) {
  AssertCompatibleTensorShape(C, R);
  Indices(C);
  C[indices] = R[indices];
}

void sisi4s::conjugate(Tensor<complex> &C) {
  Indices(C);
  CTF::Transform<complex>(std::function<void(complex &)>(
      [](complex &c) { c = std::conj(c); }))(C[indices]);
}

void sisi4s::conjugate(Tensor<double> &C) {
  // ;-)
}
