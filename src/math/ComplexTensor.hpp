#ifndef COMPLEX_TENSOR_DEFINED
#define COMPLEX_TENSOR_DEFINED

#include <math/Complex.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
/**
 * \brief Decomposes the tensor of complex elements into
 * two tensors containing the real and imaginary parts,
 * respectively.
 */
void fromComplexTensor(Tensor<complex> &c,
                       Tensor<double> &r,
                       Tensor<double> &i);

void fromComplexTensor(Tensor<double> &c, Tensor<double> &r, Tensor<double> &i);

/**
 * \brief Discards the real part of a complex tensor.
 */
void fromComplexTensor(Tensor<complex> &c, Tensor<double> &r);

/**
 * \brief Composes a tensor of complex elements
 * containing of the given tensors of real and imaginary parts.
 * Note that in this overload the imaginary part may be redistributed
 * during reading.
 */
void toComplexTensor(Tensor<double> &r, Tensor<double> &i, Tensor<complex> &c);

void toComplexTensor(Tensor<double> &r, Tensor<complex> &c);

void toComplexTensor(Tensor<double> &r, Tensor<double> &c);

void conjugate(Tensor<double> &c);

void conjugate(Tensor<complex> &c);

void conjugate(Tensor<double> &c);
} // namespace sisi4s

#endif
