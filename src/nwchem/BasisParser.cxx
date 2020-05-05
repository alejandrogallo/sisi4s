#include <nwchem/BasisParser.hpp>
#include <numeric>

using namespace pars;
using namespace nwchem;
using namespace am;

std::vector<am::AngularMomentum> am::all() {
  return { S, P, D, F, G, H, I, K };
}
size_t am::toInt (const am::AngularMomentum &a) { return a; }

am::AngularMomentum am::fromString(const std::string &a) {
  if (a == "S") return S; if (a == "P") return P;
  if (a == "D") return D; if (a == "F") return F;
  if (a == "G") return G; if (a == "H") return H;
  if (a == "I") return I; if (a == "K") return K;
  throw "I don't anderstand symbol: " + a;
}

const Regex spherical = oneOf({ "SPHERICAL", "spherical" })
          , no_line   = bof + blank + anyOf + eof
          , comment   = bof + "#" + print + anyOf
          , atom      = upper + lower + optional
          , sep       = blank + oneOrMore
          , shell_symbol  = oneOf({"S", "P", "D", "F", "G", "H", "I", "K"})
          , shell_header  = bof + capture(atom.s) + sep.s
                          + capture(shell_symbol.s)
                          + blank + anyOf
          , shell_content = bof + blank + anyOf
                          + capture(realNumber) + sep.s
                          + capture(realNumber)
          , basis_token   = oneOf({ "basis", "BASIS" })
          , basis_name    = "\"" + capture(print + oneOrMore) + "\""
          , basis_header  = group(basis_token.s) + sep.s
                          + basis_name.s + sep.s
                          + group(spherical.s)
          , basis_end     = bof + group(oneOf({ "end", "END" }))
          ;

Basis BasisSetParser::parseBasis(std::fstream &f, const std::string name) {
  std::string line;
  std::vector<Shell> shells;
  std::vector<double> coe, exp;
  std::string atom, shSymbol;

  auto addShellToBasis = [&]() {
    shells.push_back( { atom, am::fromString(shSymbol), { coe, exp } } );
    coe.resize(0); exp.resize(0);
  };

  while (std::getline(f, line)) {

    if (matches(line, shell_header)) {
      if (coe.size() and exp.size()) { addShellToBasis(); }
      atom = std::string(match[1]);
      shSymbol = std::string(match[2]);
      continue;
    }

    if (matches(line, basis_end)) {
      addShellToBasis();
      return { atom, name, shells };
    }

    if (! matches(line, shell_content) ) throw "in: " + line;

    coe.push_back(std::atof(std::string(match[2]).c_str()));
    exp.push_back(std::atof(std::string(match[1]).c_str()));

  }

  throw "Not able to parse basis, fix your input";

}

bool BasisSetParser::matches(const std::string &t, const Regex &r) {
  std::regex_match(t, match, r.r);
  return match.size() > 0;
}

BasisSet BasisSetParser::parseFile(const std::string &fileName) {
  std::fstream f(fileName);
  std::string line;
  BasisSet bs;
  if (!f.is_open()) throw "File IO error: " + fileName;

  while (std::getline(f, line)) {

    // ignore empty or comment lines
    if ( matches(line, no_line) || matches(line, comment) ) continue;

    // Parse basis whenever we find basis...
    if (matches(line, basis_header)) {
                              // v-- name
      bs.push_back(parseBasis(f, match[1].str()));
      continue;
    }

    throw "I don't understand: " + line;

  }

  return bs;

}
