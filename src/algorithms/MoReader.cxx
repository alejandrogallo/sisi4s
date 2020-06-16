#include <string>
#include <vector>
#include <algorithm>
#include <algorithms/MoReader.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>
#include <util/Parsing.hpp>
#include <util/AngularMomentum.hpp>
#include <util/BasisSet.hpp>
#include <util/Emitter.hpp>
#include <util/XyzParser.hpp>
#include <nwchem/BasisParser.hpp>
#include <nwchem/MovecsParser.hpp>
#include <turbomole/MosParser.hpp>
#include <regex>
#include <iterator>
#define LOGGER(_l) LOG(_l, "MoReader")
#define IF_GIVEN(_l, ...) if (isArgumentGiven(_l)) { __VA_ARGS__ }

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(MoReader);


std::vector< std::vector<size_t> > findShellIndices
  ( const BasisSet &bs
  , const std::vector<std::string> &atoms // atom names in our structure
  , const am::AngularMomentum &am
  ) {
  std::vector< std::vector<size_t> > indices;
  size_t atomBlock(0);

  for (const auto &a: atoms) // e.g. {N , N }
  for (const auto &b: bs) { if ( b.atom != a ) continue;

      // Basis b of atom a in the structure
      for (const auto &shell: b.shells) {
        size_t nbf(am::toInt(shell.am));
        if (shell.am != am) { // e.g. we're looking for D and get S shell
          atomBlock += nbf;
          continue;
        }
        std::vector<size_t> _indices(nbf);
        std::iota( _indices.begin()
                 , _indices.end()
                 , atomBlock
                 );
        indices.push_back(_indices);
        atomBlock += nbf;
      }

    }

  return indices;
}

struct ShellParser {
  const Regex shell_symbol = pars::oneOf({"S", "P", "D", "F", "G", "H", "I", "K"})
            , atom = pars::upper + pars::lower + pars::optional
            , sep = pars::blank + pars::oneOrMore
            , shell = pars::blank + pars::anyOf
                    // number of shells
                    + pars::capture(pars::digit + pars::oneOrMore)
                    // shell sybmol
                    + pars::capture(shell_symbol.s)
                    + pars::blank + pars::anyOf
            , basis_line = pars::capture(atom.s) + sep.s
                         + pars::capture(pars::group(shell.s)+ pars::oneOrMore)
            ;
  BasisSet parseString(const std::string &s) {
    std::smatch basisMatch, shellMatch;
    std::string S(s);
    BasisSet bs;
    while (std::regex_search(S, basisMatch, basis_line.r)) {
      std::string atomSymbol = basisMatch[1].str()
                , shellLine = basisMatch[2].str()
                ;
      std::cout << basisMatch[0].str() << std::endl;
      std::cout << "shellLine: " << shellLine << std::endl;
      std::vector<Shell> shells;
      while (std::regex_search(shellLine, shellMatch, shell.r)) {
        size_t nShells(std::atoi(shellMatch[1].str().c_str()));
        std::cout << "nshells " << nShells << std::endl;
        am::AngularMomentum am(am::fromString(shellMatch[2].str()));
        for (size_t i(0); i < nShells; i++) {
          shells.push_back({ atomSymbol, am, {} });
        }
        shellLine = shellMatch.suffix();
      }
      S = basisMatch.suffix();
      bs.push_back({atomSymbol, "Unnamed", shells});
    }
    return bs;
  }
};

std::map<std::string, std::map<std::string, std::string>>
MoReader::DEFAULT_SCALINGS =
  { { "nwchem", { {"DScaling", "1,1,1,-1,1"                  }
                , {"FScaling", "1,1,1,1,-1,1,-1"             }
                , {"GScaling", "1,1,1,1,1,-1,1,-1,1"         }
                , {"HScaling", "1,1,1,1,1,1,-1,1,-1,1,-1"    }
                , {"IScaling", "1,1,1,1,1,1,1,-1,1,-1,1,-1,1"}
                }
    }
  , { "psi4", {} }
  , { "turbomole", {} }
  };

std::map<std::string, std::map<std::string, std::string>>
MoReader::DEFAULT_REORDER =
  { { "nwchem", {} }
  , { "psi4", {} }
  , { "turbomole", { { "DReorder", "4,2,0,1,3" } } }
  };

std::vector<std::string>
MoReader::BACKENDS = {"nwchem", "psi4", "turbomole"};


void MoReader::run() {
  std::vector<std::string> args;
  const std::string fileName(getTextArgument("file"))
                  , xyz(getTextArgument("xyzStructureFile", ""))
                  , basisFile(getTextArgument("basisFile", ""))
                  , shellsFile(getTextArgument("shellsFile", ""))
                  , backend(getTextArgument("backend"))
                  ;
  const int No(getIntegerArgument("No"));

  if (std::find(BACKENDS.begin(), BACKENDS.end(), backend) == BACKENDS.end()) {
    throw "Incorrect backend value: " + backend;
  }

  // make <shell>Scaling
  const auto _s
    = [&](std::string s) { return
        std::make_pair< am::AngularMomentum
                      , std::vector<double>
                      > ( am::fromString(s)
                        , pars::parseVector<double>(
                            this->getTextArgument
                              ( s + "Scaling"
                              , this->DEFAULT_SCALINGS[backend][s + "Scaling"]
                              ))
                        );
        };
  // make <shell>Reorder
  const auto _r
    = [&](std::string s) { return
        std::make_pair< am::AngularMomentum
                      , std::vector<int>
                      > ( am::fromString(s)
                        , pars::parseVector<int>(
                            this->getTextArgument
                              ( s+ "Reorder"
                              , this->DEFAULT_REORDER[backend][s + "Reorder"]
                              ))
                        );
        };
  struct {
    const std::map<am::AngularMomentum, std::vector<double>> scaling;
    const std::map<am::AngularMomentum, std::vector<int>   > reorder;
  } transformation = { { _s("S"), _s("P"), _s("D"), _s("F")
                       , _s("G"), _s("H"), _s("I"), _s("K") }
                     , { _r("S"), _r("P"), _r("D"), _r("F")
                       , _r("G"), _r("H"), _r("I"), _r("K") }
                     };

  LOGGER(0) << "NOTE: it only works now for restricted references" << std::endl;
  LOGGER(0) << "file: " << fileName << std::endl;
  LOGGER(0) << "No: " << No << std::endl;

  BasisSet basis;
  if (basisFile.size()) {
    basis = nwchem::BasisSetParser().parseFile(basisFile);
    LOGGER(0) << "#basis in" << basisFile << ": " << basis.size() << std::endl;
  }
  if (shellsFile.size()) {
    LOGGER(0) << "basis in '" << shellsFile << "'" << std::endl;
    std::ifstream f(shellsFile);
    const std::string contents(std::istreambuf_iterator<char>(f), {});
    basis = ShellParser().parseString(contents);
  }
  for (const auto &b: basis) {
    LOGGER(0) << ":: " << b.atom << " :name " << b.name << " "
              << ":#shells" << b.shells.size() << std::endl;
  }

  std::vector<std::string> atoms;
  if (xyz.size()) {
    const auto structure(pars::XyzParser().parseFile(xyz));
    for (const auto &a: structure) atoms.push_back(a.symbol);
    LOGGER(0) << "xyzStructureFile: " << xyz << std::endl;
    LOGGER(0) << "#atoms: " << structure.size() << std::endl;
    for (const auto& a: structure) {
      LOGGER(0) << a.symbol     << ": "
                << a.position.x << ", "
                << a.position.y << ", "
                << a.position.z << std::endl;
    }
  }

  std::vector<double> mos, eigenvalues, occupations;
  size_t Np;

  if (backend == NWCHEM) {
    const nwchem::MovecReader movec(fileName);
    mos = movec.mos;
    eigenvalues = movec.eigenvalues;
    occupations = movec.occupations;
    Np = movec.Np;
  } else if (backend == TURBOMOLE) {
    const tmole::MosParser movec(fileName);
    mos = movec.mos;
    eigenvalues = movec.eigenvalues;
    occupations = movec.occupations;
    Np = movec.Np;
  }

  for (const auto &p: transformation.scaling) {
    const auto am = am::toInt(p.first);
    const auto& scaling = p.second;
    if (scaling.size()) {
      LOGGER(0) << "doing scaling of am: "  << am << std::endl;
      LOGGER(0) << "You should have "  << am << " numbers" << std::endl;
      assert(scaling.size() == am);
      auto indicesVector(findShellIndices(basis, atoms, p.first));
      // go through every state (column of mos)
      for (size_t j(0); j<Np; j++) {
      for (auto const& indices: indicesVector) {
        assert(indices.size() == am);
      for (size_t ii(0); ii<am; ii++) {
        assert(indices[ii] < Np);
        mos[ indices[ii] + j * Np ] *= scaling[ii];
      }
      }
      }
      LOGGER(0) << "done with "  << am << std::endl;
    }
  }

  for (const auto &p: transformation.reorder) {
    const auto am = am::toInt(p.first);
    const auto& reorder = p.second;
    if (reorder.size()) {
      LOGGER(0) << "doing reorder of am: "  << am << std::endl;
      LOGGER(0) << "You should have "  << am << " numbers" << std::endl;
      assert(reorder.size() == am);
      auto indicesVector(findShellIndices(basis, atoms, p.first));
      // go through every state (column of mos)
      for (size_t j(0); j<Np; j++) {
      for (auto const& indices: indicesVector) {
        assert(indices.size() == am);
        // save a backup of the indices thunk to be settable
        std::vector<double> backup(indices.size());
        for (size_t ii(0); ii<am; ii++) backup[ii] = mos[indices[ii] + j*Np];
      for (size_t ii(0); ii<am; ii++) {
        assert(indices[ii] < Np);
        mos[ indices[ii] + j * Np ] = backup[reorder[ii]];
      }
      }
      }
      LOGGER(0) << "done with "  << am << std::endl;
    }
  }

  std::vector<int> pp(2, Np), syms(2, NS), o(1, No), v(1, Np - No);
  std::vector<int64_t> ids;
  const int rank_m = int(Cc4s::world->rank == 0); // rank mask


  IF_GIVEN("HoleEigenEnergies",
    ids.resize(rank_m * o[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsi(new CTF::Tensor<double>(1, o.data(), syms.data(), *Cc4s::world));
    epsi->write(ids.size(), ids.data(), eigenvalues.data());
    allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  )

  IF_GIVEN("ParticleEigenEnergies",
    ids.resize(rank_m * v[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsa(new CTF::Tensor<double>(1, v.data(), syms.data(), *Cc4s::world));
    epsa->write(ids.size(), ids.data(), eigenvalues.data() + No);
    allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);
  )

  IF_GIVEN("OccupationNumbers",
    ids.resize(rank_m * Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto os(new CTF::Tensor<double>(1, pp.data(), syms.data(), *Cc4s::world));
    os->write(ids.size(), ids.data(), occupations.data());
    allocatedTensorArgument<double>("OccupationNumbers", os);
  )

  if (isArgumentGiven("OrbitalCoefficients")) {
    ids.resize(rank_m * Np*Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto coef(new CTF::Tensor<double>(2, pp.data(), syms.data(), *Cc4s::world));
    coef->write(ids.size(), ids.data(), mos.data());
    allocatedTensorArgument<double>("OrbitalCoefficients", coef);
  }

}
