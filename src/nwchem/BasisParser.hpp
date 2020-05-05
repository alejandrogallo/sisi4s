#ifndef _NWCHEM_BASIS_PARSER_
#define _NWCHEM_BASIS_PARSER_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <util/Parsing.hpp>
#include <regex>


namespace nwchem {

  using namespace pars;

  struct ContractedGaussian {
    const std::vector<double> coefficients, exponents;
    size_t size() const { return coefficients.size(); }
  };

  namespace am {
    enum AngularMomentum { S = 1, P = 3 , D = 5 , F = 7
                         , G = 9, H = 11, I = 13, K = 15
                         };
    std::vector<AngularMomentum> all();
    size_t toInt (const AngularMomentum &);
    AngularMomentum fromString(const std::string &);
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
    bool matches(const std::string &t, const Regex &r);
    Basis parseBasis (std::fstream &f, const std::string name);
    BasisSet parseFile(const std::string &fileName);
  };

}

#endif
