#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include <util/CTF.hpp>

namespace sisi4s {

  template <typename F>
  using Tensor = CTF::Tensor<F>;

  template <typename F>
  using Matrix = CTF::Matrix<F>;

}  // namespace sisi4s
#endif
