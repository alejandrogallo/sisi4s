#ifndef COULOMBINTEGRALS_HPP_
#define COULOMBINTEGRALS_HPP_

#include <map>

#include <util/Tensor.hpp>

#define _DEFINE_SETTER(__type, __name, no, nv)                                 \
  CoulombIntegrals<F> &with_##__name(__type t) {                               \
    if (no > 0) { No = t->lens[no]; }                                          \
    if (nv > 0) { Nv = t->lens[nv]; }                                          \
    __name##_ptr = t;                                                          \
    return *this;                                                              \
  }                                                                            \
  __type __name() {                                                            \
    if (!__name##_ptr == nullptr) return __name##_ptr;                         \
    throw "TODO: Automatic calculation from Γ not yet implemented!";           \
  }                                                                            \
  __type __name##_ptr = nullptr

namespace sisi4s {

template <typename F>
class CoulombIntegrals {
public:
  using tensor_type = Tensor<F> *;

  size_t No = 0, Nv = 0;
  Tensor<sisi4s::complex> *gamma = nullptr;

  CoulombIntegrals(size_t no, size_t nv, Tensor<sisi4s::complex> *g)
      : No(no)
      , Nv(nv)
      , gamma(g) {}

  CoulombIntegrals() {}

  _DEFINE_SETTER(tensor_type, hhhh, 0, -1);
  _DEFINE_SETTER(tensor_type, hhhp, 0, 3);
  _DEFINE_SETTER(tensor_type, hhph, 0, 2);
  _DEFINE_SETTER(tensor_type, hhpp, 0, 2);
  _DEFINE_SETTER(tensor_type, hphh, 0, 1);
  _DEFINE_SETTER(tensor_type, hphp, 0, 1);
  _DEFINE_SETTER(tensor_type, hpph, 0, 1);
  _DEFINE_SETTER(tensor_type, hppp, 0, 1);
  _DEFINE_SETTER(tensor_type, phhh, 1, 0);
  _DEFINE_SETTER(tensor_type, phhp, 1, 0);
  _DEFINE_SETTER(tensor_type, phph, 1, 0);
  _DEFINE_SETTER(tensor_type, phpp, 1, 0);
  _DEFINE_SETTER(tensor_type, pphh, 2, 0);
  _DEFINE_SETTER(tensor_type, pphp, 2, 0);
  _DEFINE_SETTER(tensor_type, ppph, 3, 0);
  _DEFINE_SETTER(tensor_type, pppp, -1, 1);
};

} // namespace sisi4s

#undef _DEFINE_SETTER
#endif
