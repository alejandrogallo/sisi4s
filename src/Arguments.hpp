#ifndef ARGUMENTS_HPP_
#define ARGUMENTS_HPP_

#include <string>
#include <vector>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <util/Format.hpp>
#include <util/Tensor.hpp>
#include <math/Complex.hpp>

#include <NewData.hpp>
namespace sisi4s {

class Arguments {
public:
  using ArgumentName = std::string;
  using Index = std::string;
  std::map<ArgumentName, Index> arguments;

  static std::string normalize_name(std::string const &name);
  void push(std::string const &name, std::string const &db_index);
  bool present(std::string const &name);
  data::StorePair get_data(std::string const &name);
  std::string get_var(std::string const &name);

  template <typename F>
  bool is_of_type(std::string const &name) {
    return data::istype<F>(get_data(normalize_name(name)).first);
  }

  template <typename F>
  F *getptr(std::string const &name) {
    const std::string nname = normalize_name(name);
    auto const pair = get_data(nname);
    const std::string db_index = pair.first;
    return data::get<F>(db_index);
  }

  template <typename F>
  F get(std::string const &name) {
    const std::string nname = normalize_name(name);
    F value = *getptr<F>(nname);
    EMIT() << YAML::Key << nname << YAML::Value << value;
    return value;
  }

  template <typename F>
  void set(std::string const &name, F const &value) {
    F *ptr = getptr<F>(normalize_name(name));
    *ptr = value;
  }

  template <typename F>
  void set_force(std::string const &name, F const &value) {
    const std::string nname = normalize_name(name);
    auto const pair = get_data(nname);
    const std::string db_index = pair.first;
    set<F>(nname, value);
    data::declare_type<F>(db_index);
  }
};

} // namespace sisi4s

#endif
