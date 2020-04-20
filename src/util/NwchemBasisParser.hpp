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

  struct Shell {
    const std::string atom;
    const std::string symbol;
    const ContractedGaussian g;
  };

  struct Basis {
    const std::string atom;
    const std::vector<Shell> shells;
  };

  typedef std::vector<Basis> BasisSet;

  struct BasisSetParser {

    std::smatch match;

    const
    std::string spherical = oneOf({ "SPHERICAL" })
              , emptyLine = bof + blank + anyOf + eof
              , shell_symbol = oneOf({"S", "P", "D", "F", "G", "H", "I", "K"})
              , comment = bof + "#" + print + anyOf
              , sep = blank + oneOrMore // any number > 1 of spaces or tabs
              , atom = upper + lower + optional
              , shell_header = bof + capture(atom) + sep
                             + capture(shell_symbol)
              , shell_content = bof + blank + anyOf
                              + capture(realNumber) + sep + capture(realNumber)
              , basis_header = "basis" + sep + print + oneOrMore
                             + sep + spherical
              , basis_end = oneOf({ "end", "END" })
              ;

    bool matches(const std::string &t, const std::string &m) {
      // creating a regex every time is not my bottleneck, so..
      std::regex_match(t, match, std::regex(m));
      return match.size() > 0;
    }

    Basis parseBasis (std::fstream &f) {
      std::string line;
      std::vector<Shell> shells;
      std::vector<double> coe, exp;
      std::string atom, shSymbol;

      auto addShellToBasis = [&]() {
          shells.push_back( { atom, shSymbol, { coe, exp } } );
          coe.resize(0); exp.resize(0);
      };

      while (std::getline(f, line)) {

        if (matches(line, shell_header)) {
          atom = std::string(match[1]);
          shSymbol = std::string(match[2]);
          if (coe.size() and exp.size()) { addShellToBasis(); }
          continue;
        }

        if (matches(line, basis_end)) {
          addShellToBasis();
          return { atom, shells };
        }

        if (! matches(line, shell_content) ) throw EXCEPTION("in: " + line);

        coe.push_back(std::atof(std::string(match[1]).c_str()));
        exp.push_back(std::atof(std::string(match[2]).c_str()));

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
        if ( matches(line, emptyLine) || matches(line, comment) ) continue;

        // We don't really care about the `basis "blasdf" SPHERICAL` header
        if (matches(line, basis_header)) {
          std::cout << line << std::endl;
          bs.push_back(parseBasis(f));
          continue;
        }

        throw EXCEPTION("I don't understand: " + line);

      }

      return bs;

    }

  };

}

#endif
