#include <util/Tensor.hpp>

namespace sisi4s {
namespace data {

#define INSTANTIATE(type)                                                      \
  std::string Namer<Tensor<type> *>::name() {                                  \
    return "Tensor of " + Namer<type>::name();                                 \
  }
INSTANTIATE(double)
INSTANTIATE(sisi4s::complex)
#undef INSTANTIATE
} // namespace data
} // namespace sisi4s
