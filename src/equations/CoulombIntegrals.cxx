#include <equations/CoulombIntegrals.hpp>

namespace sisi4s {

// TODO
// template <typename F>
// static typename CoulombIntegrals<F>::tensor_type
// make_integral(CoulombIntegrals<F> &vpqrs, std::string const &name) {
//   auto map = vpqrs.integral_map;
//   if (map[name] != nullptr) return map[name];
//   return nullptr;
// }

// template <typename F>
// typename CoulombIntegrals<F>::tensor_type CoulombIntegrals<F>::hhhh() {
//   return make_integral(*this, "hhhh");
// }

template <>
class CoulombIntegrals<double>;
template <>
class CoulombIntegrals<complex>;

} // namespace sisi4s
