#ifndef _NWCHEM_BASIS_PARSER_
#define _NWCHEM_BASIS_PARSER_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <util/Parsing.hpp>
#include <util/AngularMomentum.hpp>
#include <util/BasisSet.hpp>
#include <regex>
#include <numeric>

namespace nwchem {

using namespace pars;

struct BasisSetParser {
  std::smatch match;
  bool matches(const std::string &t, const Regex &r);
  Basis parseBasis(std::fstream &f, const std::string name);
  BasisSet parseFile(const std::string &fileName);
};

} // namespace nwchem

#endif
