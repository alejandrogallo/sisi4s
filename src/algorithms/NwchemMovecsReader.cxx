#include <string>
#include <vector>
#include <algorithm>
#include <algorithms/NwchemMovecsReader.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>
#include <util/Parsing.hpp>
#include <util/Emitter.hpp>
#include <util/XyzParser.hpp>
#include <nwchem/BasisParser.hpp>
#include <nwchem/MovecsParser.hpp>
#include <regex>
#define LOGGER(_l) LOG(_l, "NwchemMovecsReader")
#define IF_GIVEN(_l, ...) if (isArgumentGiven(_l)) { __VA_ARGS__ }

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(NwchemMovecsReader);


std::vector< std::vector<size_t> > findShellIndices
  ( const nwchem::BasisSet &bs
  , const std::vector<std::string> &atoms // atom names in our structure
  , const nwchem::am::AngularMomentum &am
  ) {
  std::vector< std::vector<size_t> > indices;
  size_t atomBlock(0);

  for (const auto &a: atoms) // e.g. {N , N }
  for (const auto &b: bs) { if ( b.atom != a ) continue;

      // Basis b of atom a in the structure
      for (const auto &shell: b.shells) {
        size_t nbf(nwchem::am::toInt(shell.am));
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

void NwchemMovecsReader::run() {
  std::vector<std::string> args;
  const std::string fileName(getTextArgument("file"))
                  , xyz(getTextArgument("xyzStructureFile", ""))
                  , basisFile(getTextArgument("basisFile", ""))
                  ;
  const int No(getIntegerArgument("No"));
  std::map<std::string, std::string>
    nwchemScalings = { {"DScaling", "1,1,1,-1,1"}
                     , {"FScaling", "1,1,1,1,-1,1,-1"}
                     , {"Gscaling", "1,1,1,1,1,-1,1,-1,1"}
                     };

  // make <shell>Scaling
  const auto _s
    = [&](std::string s) { return
        std::make_pair< nwchem::am::AngularMomentum
                      , std::vector<double>
                      > ( nwchem::am::fromString(s)
                        , pars::parseVector<double>(
                            this->getTextArgument
                              ( s + "Scaling"
                              , nwchemScalings[s + "Scaling"]
                              ))
                        );
        };
  // make <shell>Reorder
  const auto _r
    = [&](std::string s) { return
        std::make_pair< nwchem::am::AngularMomentum
                      , std::vector<int>
                      > ( nwchem::am::fromString(s)
                        , pars::parseVector<int>(
                            this->getTextArgument(s+ "Reorder", ""))
                        );
        };
  struct {
    const std::map<nwchem::am::AngularMomentum, std::vector<double>> scaling;
    const std::map<nwchem::am::AngularMomentum, std::vector<int>   > reorder;
  } transformation = { { _s("S"), _s("P"), _s("D"), _s("F")
                       , _s("G"), _s("H"), _s("I"), _s("K") }
                     , { _r("S"), _r("P"), _r("D"), _r("F")
                       , _r("G"), _r("H"), _r("I"), _r("K") }
                     };

  LOGGER(0) << "NOTE: it only works now for restricted references" << std::endl;
  LOGGER(0) << "file: " << fileName << std::endl;
  LOGGER(0) << "No: " << No << std::endl;

  nwchem::BasisSet basis;
  if (basisFile.size()) {
    basis = nwchem::BasisSetParser().parseFile(basisFile);
    LOGGER(0) << "#basis in"
              << basisFile << ": " << basis.size() << std::endl;
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

  const nwchem::MovecReader movec(fileName);
  std::vector<double> mos(movec.Np * movec.Np);
  mos = movec.mos;

  for (const auto &p: transformation.scaling) {
    const auto am = nwchem::am::toInt(p.first);
    const auto& scaling = p.second;
    const auto& Np = movec.Np;
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
    const auto am = nwchem::am::toInt(p.first);
    const auto& reorder = p.second;
    const auto& Np = movec.Np;
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

  std::vector<int> pp(2, movec.Np), syms(2, NS), o(1, No), v(1, movec.Np - No);
  std::vector<int64_t> ids;
  const int rank_m = int(Cc4s::world->rank == 0); // rank mask


  IF_GIVEN("HoleEigenEnergies",
    ids.resize(rank_m * o[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsi(new CTF::Tensor<double>(1, o.data(), syms.data(), *Cc4s::world));
    epsi->write(ids.size(), ids.data(), movec.eigenvalues.data());
    allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  )

  IF_GIVEN("ParticleEigenEnergies",
    ids.resize(rank_m * v[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsa(new CTF::Tensor<double>(1, v.data(), syms.data(), *Cc4s::world));
    epsa->write(ids.size(), ids.data(), movec.eigenvalues.data() + No);
    allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);
  )

  IF_GIVEN("OccupationNumbers",
    ids.resize(rank_m * movec.Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto os(new CTF::Tensor<double>(1, pp.data(), syms.data(), *Cc4s::world));
    os->write(ids.size(), ids.data(), movec.occupations.data());
    allocatedTensorArgument<double>("OccupationNumbers", os);
  )

  if (isArgumentGiven("OrbitalCoefficients")) {
    ids.resize(rank_m * movec.Np*movec.Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto coef(new CTF::Tensor<double>(2, pp.data(), syms.data(), *Cc4s::world));
    coef->write(ids.size(), ids.data(), mos.data());
    allocatedTensorArgument<double>("OrbitalCoefficients", coef);
  }

}
