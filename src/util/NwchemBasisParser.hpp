#ifndef _NWCHEM_BASIS_PARSER_
#define _NWCHEM_BASIS_PARSER_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <util/Parsing.hpp>
#include <regex>

using namespace pars;


namespace nwchem {

  struct ContractedGaussian {
    const std::vector<double> coefficients, exponents;
  };

  namespace am {
    enum AngularMomentum { S = 1, P = 3 , D = 5 , F = 7
                         , G = 9, H = 11, I = 13, K = 15
                         };
    size_t toInt (const AngularMomentum &am) { return am; }
    AngularMomentum fromString(const std::string &am) {
      if (am == "S") return S; if (am == "P") return P;
      if (am == "D") return D; if (am == "F") return F;
      if (am == "G") return G; if (am == "H") return H;
      if (am == "I") return I; if (am == "K") return K;
      throw EXCEPTION("I don't anderstand symbol: " + am);
    }
  }

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
      return
      std::accumulate( shells.begin()
                     , shells.end()
                     , 0
                     , [&](size_t i, const Shell &s){ return i + s.size(); }
                     );
    }
  };

  typedef std::vector<Basis> BasisSet;

  struct BasisSetParser {

    std::smatch match;

    const Regex spherical = oneOf({ "SPHERICAL", "spherical" })
              , no_line   = bof + blank + anyOf + eof
              , comment   = bof + "#" + print + anyOf
              , atom      = upper + lower + optional
              , sep       = blank + oneOrMore
              , shell_symbol  = oneOf({"S", "P", "D", "F", "G", "H", "I", "K"})
              , shell_header  = bof + capture(atom.s) + sep.s
                              + capture(shell_symbol.s)
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

    bool matches(const std::string &t, const Regex &r) {
      std::regex_match(t, match, r.r);
      return match.size() > 0;
    }

    Basis parseBasis (std::fstream &f, const std::string name) {
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

        if (! matches(line, shell_content) ) throw EXCEPTION("in: " + line);

        coe.push_back(std::atof(std::string(match[2]).c_str()));
        exp.push_back(std::atof(std::string(match[1]).c_str()));

      }

      throw EXCEPTION("Not able to parse basis, fix your input");

    }

    BasisSet parseFile(const std::string &fileName) {
      std::fstream f(fileName);
      std::string line;
      BasisSet bs;
      if (!f.is_open()) throw EXCEPTION("File IO error: " + fileName);

      while (std::getline(f, line)) {

        // ignore empty or comment lines
        if ( matches(line, no_line) || matches(line, comment) ) continue;

        // Parse basis whenever we find basis...
        if (matches(line, basis_header)) {
                                  // v-- name
          bs.push_back(parseBasis(f, match[1].str()));
          continue;
        }

        throw EXCEPTION("I don't understand: " + line);

      }

      return bs;

    }

  };

}

#endif
