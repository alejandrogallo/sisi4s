#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include <util/CTF.hpp>
#include <math/Complex.hpp>
#include <NewData.hpp>

namespace sisi4s {

template <typename F>
using Tensor = CTF::Tensor<F>;

using TensorIndex = int64_t;

template <typename F>
using Matrix = CTF::Matrix<F>;

namespace data {

template <typename F>
class Namer<Tensor<F> *> {
public:
  static std::string name() { return "Tensor of " + Namer<F>::name(); }
};

#define INSTANTIATE(type)                                                      \
  template <>                                                                  \
  class Namer<Tensor<type> *> {                                                \
  public:                                                                      \
    static std::string name();                                                 \
  };
INSTANTIATE(double)
INSTANTIATE(sisi4s::complex)
#undef INSTANTIATE

} // namespace data

} // namespace sisi4s
#endif
