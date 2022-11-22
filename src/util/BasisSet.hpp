#ifndef ___BASIS_SET_DEFINITIONS_
#define ___BASIS_SET_DEFINITIONS_

#include <vector>
#include <util/AngularMomentum.hpp>
#include <numeric>

struct ContractedGaussian {
  const std::vector<double> coefficients, exponents;
  size_t size() const { return coefficients.size(); }
};

struct Shell {
  const std::string atom;
  const am::AngularMomentum am;
  const ContractedGaussian g;
  size_t size() const { return g.size(); }
};

struct Basis {
  const std::string atom;
  const std::string name;
  const std::vector<Shell> shells;
  size_t size() const { return shells.size(); }
  size_t nbf() const {
    return std::accumulate(
        shells.begin(),
        shells.end(),
        0,
        [&](size_t i, const Shell &s) { return i + s.size(); });
  }
};

typedef std::vector<Basis> BasisSet;

#endif
