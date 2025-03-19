#ifndef SPECCHAIN_HPP
#define SPECCHAIN_HPP

#include <AlgorithmInputSpec.hpp>

namespace sisi4s {
namespace spec {

class SpecChain : public Spec {
public:
  const std::vector<Spec *> subspecs;
  Spec *valid_spec = nullptr;
  SpecChain(const std::vector<Spec *> _subspecs)
      : subspecs(_subspecs) {}
  virtual std::vector<std::string> autodoc() override {
    std::vector<std::string> r;
    for (auto const &s : subspecs) {
      for (auto const &i : s->autodoc()) { r.push_back(i); }
    }
    return r;
  }
  virtual std::string type_user_name() override {
    std::stringstream out;
    for (auto const &s : subspecs) { out << s->type_user_name() << " | "; }
    return out.str();
  }
  virtual void parse(std::string const &val) override {
    for (auto const &s : subspecs) { s->parse(val); }
  }
  virtual std::string commit() override {
    validate();
    return valid_spec->commit();
  }
  virtual bool validate() override {
    for (auto const &s : subspecs) {
      if (s->validate()) {
        valid_spec = s;
        return true;
      }
    }
    return false;
  }
  virtual std::vector<std::string> warnings(std::string const &v) {
    IGNORABLE(v);
    return {"None of the specs match the given input"};
  }
};

#define SPEC_CHAIN(doc, ...)                                                   \
  ((new sisi4s::spec::SpecChain(                                               \
        std::vector<sisi4s::spec::Spec *>{__VA_ARGS__}))                       \
       ->with_doc(doc)                                                         \
       ->with_default(false))

} // namespace spec
} // namespace sisi4s

#endif
