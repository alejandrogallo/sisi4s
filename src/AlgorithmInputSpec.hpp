#ifndef ALGORITHMINPUTSPEC_HPP_
#define ALGORITHMINPUTSPEC_HPP_

#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <functional>
#include <iostream>

#include <NewData.hpp>
#include <util/Ignorable.hpp>

#include <yaml-cpp/yaml.h>

namespace sisi4s {

namespace spec {

class Spec {
public:
  // todo: required?
  // todo: warnings
  virtual std::vector<std::string> warnings(std::string const &v) = 0;
  virtual std::string commit() = 0;
  virtual void parse(std::string const &) = 0;
  virtual bool validate() = 0;
  std::string doc;
  Spec *with_doc(std::string const &d) {
    doc = d;
    return this;
  }
  bool has_default = false;
  Spec *with_default(bool v) {
    has_default = v;
    return this;
  }
  bool required = false;
  virtual Spec *require() {
    required = true;
    return this;
  }
};

template <typename F>
class YAMLSpec : public Spec {
public:
  F value;
  virtual std::string commit() override { return data::put<F>(value).first; }
  virtual void parse(std::string const &v) override {
    YAML::Node n;
    n["v"] = v;
    try {
      value = n["v"].as<F>();
    } catch (YAML::TypedBadConversion<F> const &c) {
      std::cout << "We could not parse the value '" << v << "'" << std::endl;
      std::cout << c.what() << std::endl;
      throw "";
    }
  }
  virtual bool validate() override { return true; }
  virtual std::vector<std::string> warnings(std::string const &v) {
    IGNORABLE(v);
    return {"A value of type F"};
  }
};

template <typename F>
class InVariable : public Spec {
public:
  bool found_variable = false;
  std::string db_index;
  virtual void parse(std::string const &var_name) override {
    const auto pair = data::find_by_name(var_name);
    found_variable = !data::null(pair);
    db_index = pair.first;
  }
  virtual bool validate() override { return found_variable; }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    if (found_variable) {
      return {};
    } else {
      return {"The Variable '" + v + "' was not previously declared"};
    }
  }
  virtual std::string commit() override { return db_index; }
};
#define SPEC_VARIN(doc, type)                                                  \
  (new sisi4s::spec::InVariable<type>())->with_doc(doc)

template <typename F>
class OutVariable : public InVariable<F> {
public:
  bool found_variable = false;
  std::string varname;
  virtual void parse(std::string const &var_name) override {
    const auto pair = data::find_by_name(var_name);
    found_variable = !data::null(pair);
    varname = var_name;
  }
  virtual bool validate() override { return !found_variable; }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    if (!found_variable) {
      return {};
    } else {
      return {"Nameclash: The Variable '" + v + "' was previously declared"};
    }
  }
  virtual std::string commit() override {
    auto pair = data::put<F>(F());
    std::string db_index = pair.first;
    auto data = data::getraw(db_index);
    data->name = varname;
    return db_index;
  }
};
#define SPEC_VAROUT(doc, type)                                                 \
  (new sisi4s::spec::OutVariable<type>())->with_doc(doc)

template <typename F>
class OneOf : public YAMLSpec<F> {
public:
  const std::vector<F> options;
  OneOf(std::vector<F> const &options_)
      : options(options_) {
    this->value = options[0];
  }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    std::stringstream s;
    s << v << " is not one of: <  ";
    for (auto const &o : options) s << o << "  ";
    s << ">";
    return {s.str()};
  }
  virtual bool validate() override {
    return std::find(options.begin(), options.end(), this->value)
        != options.end();
  }
};
#define SPEC_ONE_OF(doc, type, ...)                                            \
  (new sisi4s::spec::OneOf<type>({__VA_ARGS__}))                               \
      ->with_doc(doc)                                                          \
      ->with_default(true)

template <typename F>
class Value : public YAMLSpec<F> {
public:
  Value() {}
  Value(const F def) { this->value = def; }
  virtual void parse(std::string const &v) override {
    if (v.size()) { YAMLSpec<F>::parse(v); }
  }
  virtual bool validate() override { return true; }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    IGNORABLE(v);
    std::stringstream s;
    s << "Shoulde be a value";
    return {s.str()};
  }
};
#define SPEC_VALUE_DEF(doc, type, ...)                                         \
  (new sisi4s::spec::Value<type>(__VA_ARGS__))                                 \
      ->with_doc(doc)                                                          \
      ->with_default(true)
#define SPEC_VALUE(doc, type)                                                  \
  (new sisi4s::spec::Value<type>())->with_doc(doc)->with_default(false)

template <typename F>
class Satisfies : public Value<F> {
public:
  using Lambda = std::function<bool(F const &)>;
  const Lambda lambda;
  const std::string description;
  Satisfies(std::string const &description_, Lambda const &lambda_)
      : lambda(lambda_)
      , description(description_) {}
  virtual bool validate() override { return lambda(this->value); }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    std::stringstream s;
    s << "Value " << v << " does not satisfy condition: " << description;
    return {s.str()};
  }
};
#define SPEC_SATISFIES(doc, type, description, lambda)                         \
  (new sisi4s::spec::Satisfies<type>(                                          \
       description,                                                            \
       sisi4s::spec::Satisfies<type>::Lambda(lambda)))                         \
      ->with_doc(doc)                                                          \
      ->with_default(false)
#define SPEC_POSITIVE(doc, type)                                               \
  SPEC_SATISFIES(doc, type, "Greater than 0", [](type const &v) -> bool {      \
    return v > 0;                                                              \
  })
#define SPEC_NEGATIVE(doc, type)                                               \
  SPEC_SATISFIES(doc, type, "Less than 0", [](type const &v) -> bool {         \
    return v < 0;                                                              \
  })

template <typename F>
class Range : public YAMLSpec<F> {
public:
  F low, high, value;
  Range(F low_, F high_)
      : low(low_)
      , high(high_) {}
  virtual bool validate() override { return value < high && low < value; }
  virtual std::vector<std::string> warnings(std::string const &v) override {
    std::stringstream s;
    s << "Should be in the interval between " << low << " and " << high;
    return {s.str()};

    std::vector<std::string> result = this->warnings(v);
    result.push_back(s.str());
    return result;
  }
};
#define SPEC_RANGE(doc, type, high, low)                                       \
  (new sisi4s::spec::Range<type>(high, low))->with_doc(doc)

} // namespace spec

class AlgorithmInputSpec {
public:
  using Spec = std::map<std::string, spec::Spec *>;
  struct IOSpec {
    Spec in, out;
  };
  using Map = std::map<std::string, IOSpec>;
  static Map map;
  AlgorithmInputSpec(std::string const &name, Spec spec_in, Spec spec_out) {
    map[name] = {spec_in, spec_out};
  }
};

} // namespace sisi4s

#define SPEC_IN(...)                                                           \
  { __VA_ARGS__ }
#define SPEC_OUT(...)                                                          \
  { __VA_ARGS__ }
#define DEFSPEC(NAME, ...) AlgorithmInputSpec NAME::spec(#NAME, __VA_ARGS__)

#endif
