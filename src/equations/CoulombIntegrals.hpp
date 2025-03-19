#ifndef COULOMBINTEGRALS_HPP_
#define COULOMBINTEGRALS_HPP_

#include <map>

#include <util/Tensor.hpp>

#define _DEFINE_SETTER(__type, __name, __default)                              \
  CoulombIntegrals<F> &with_##__name(__type t) {                               \
    __name##_ptr = t;                                                          \
    integral_map[#__name] = t;                                                 \
    return *this;                                                              \
  }                                                                            \
  __type __name() {                                                            \
    if (!__name##_ptr == nullptr) return __name##_ptr;                         \
    throw "TODO: Automatic calculation from Γ not yet implemented!";           \
  }                                                                            \
  __type __name##_ptr = __default

namespace sisi4s {

template <typename F>
class CoulombIntegrals {
public:
  using tensor_type = Tensor<F> *;
  std::map<std::string, tensor_type> integral_map;

  size_t No, Nv;
  Tensor<sisi4s::complex> *gamma;

  CoulombIntegrals(size_t no, size_t nv, Tensor<sisi4s::complex> *g)
      : No(no)
      , Nv(nv)
      , gamma(g) {}

  CoulombIntegrals() {}

  _DEFINE_SETTER(tensor_type, hhhh, nullptr);
  _DEFINE_SETTER(tensor_type, hhhp, nullptr);
  _DEFINE_SETTER(tensor_type, hhph, nullptr);
  _DEFINE_SETTER(tensor_type, hhpp, nullptr);
  _DEFINE_SETTER(tensor_type, hphh, nullptr);
  _DEFINE_SETTER(tensor_type, hphp, nullptr);
  _DEFINE_SETTER(tensor_type, hpph, nullptr);
  _DEFINE_SETTER(tensor_type, hppp, nullptr);
  _DEFINE_SETTER(tensor_type, phhh, nullptr);
  _DEFINE_SETTER(tensor_type, phhp, nullptr);
  _DEFINE_SETTER(tensor_type, phph, nullptr);
  _DEFINE_SETTER(tensor_type, phpp, nullptr);
  _DEFINE_SETTER(tensor_type, pphh, nullptr);
  _DEFINE_SETTER(tensor_type, pphp, nullptr);
  _DEFINE_SETTER(tensor_type, ppph, nullptr);
  _DEFINE_SETTER(tensor_type, pppp, nullptr);
};

} // namespace sisi4s

#undef _DEFINE_SETTER
#endif
