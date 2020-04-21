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
#include <util/Emitter.hpp>
#include <util/Parsing.hpp>
#include <util/XyzParser.hpp>
#include <util/NwchemBasisParser.hpp>
#include <regex>
#define LOGGER(_l) LOG(_l, "NwchemMovecsReader")
#define IF_GIVEN(_l, ...) if (isArgumentGiven(_l)) { __VA_ARGS__ }

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(NwchemMovecsReader);

template <typename F>
std::vector<F> readFortranChunk(std::fstream &file) {
  uint32_t fortran_section, copy;
  std::vector<F> buffer;

  file.read((char*)&fortran_section, sizeof(int32_t));
    copy = fortran_section;
    buffer.resize(fortran_section / sizeof(F));
    file.read((char*)buffer.data(), fortran_section);
  file.read((char*)&fortran_section, sizeof(int32_t));
  assert(copy == fortran_section);

  return buffer;
}


struct MovecReader {
  size_t nsets, Np;
  std::vector<size_t> nmo;
  std::vector<double> occupations, eigenvalues, mos;
  size_t iset = 0; // This allows to skip the first set, dont need this yet
  MovecReader(const std::string &fileName) {
    std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    readFortranChunk<char>(file); // convergence info
    readFortranChunk<char>(file); // scftype
    readFortranChunk<char>(file); // lentit
    LOGGER(0) << "title: " << readFortranChunk<char>(file).data() << std::endl;
    readFortranChunk<char>(file).data(); // lenbas
    LOGGER(0) << "basis: " << readFortranChunk<char>(file).data() << std::endl;

    nsets = *readFortranChunk<int64_t>(file).data(); // nsets
    Np = *readFortranChunk<int64_t>(file).data();   // Np

    LOGGER(0) << "nsets: " << nsets << std::endl;
    LOGGER(0) << "Np: " << Np << std::endl;

    // read nmos
    for (size_t i=0; i<nsets; i++) {
      int64_t buf = *readFortranChunk<int64_t>(file).data();
      nmo.push_back(buf);
      LOGGER(0) << "nmo[" << i << "]: " << nmo[0] << std::endl;
    }

    // loop which reject MOs from the set
    for (size_t i = 0; i < iset;  i++) {
      readFortranChunk<char>(file);
      readFortranChunk<char>(file);
      for (uint64_t j=0; j<nmo[i]; j++) { //do i = 1, nmo(jset) read(unitno)
        readFortranChunk<char>(file);
      }
    }

    for (size_t s(0); s < nsets; s++) { // do s = 1, nsets

      // Read occupation numbers
      occupations = readFortranChunk<double>(file);
      LOGGER(0) << "#occupations read: " << occupations.size() << std::endl;
      assert(occupations.size() == Np);

      // Eigenvaluess [read(unitno) (evals(j), j=1,Np)]
      eigenvalues = readFortranChunk<double>(file);
      LOGGER(0) << "#eigenvalues read: " <<  eigenvalues.size() << std::endl;
      assert(eigenvalues.size() == Np);

      mos.resize(Np * Np);
      for (size_t i(0); i < nmo[iset]; i++) { // do i = 1, nmo(iset)
        std::vector<double> mos_buff(Np);
        mos_buff = readFortranChunk<double>(file);
        assert(mos_buff.size() == Np);
        for (size_t j(0); j < mos_buff.size(); j++)
          mos[j + Np*i] = mos_buff[j];
      }
    }
  }
};

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
                  , xyz(getTextArgument("xyzStructureFile"))
                  , basisFile(getTextArgument("basisFile"))
                  ;
  const int No(getIntegerArgument("No"));

  // get <shell>Scaling
  const auto _s
    = [&](std::string s) { return
        std::make_pair< nwchem::am::AngularMomentum
                      , std::vector<double>
                      > ( nwchem::am::fromString(s)
                        , parseVector<double>(
                            this->getTextArgument(s+ "Scaling", ""))
                        );
        };
  // get <shell>Reorder
  const auto _r
    = [&](std::string s) { return
        std::make_pair< nwchem::am::AngularMomentum
                      , std::vector<int>
                      > ( nwchem::am::fromString(s)
                        , parseVector<int>(
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

  const auto basis(nwchem::BasisSetParser().parseFile(basisFile));
  LOGGER(0) << "#basis in " << basisFile << ": " << basis.size() << std::endl;

  const auto structure(XyzParser().parseFile(xyz));
  std::vector<std::string> atoms;
  for (const auto &a: structure) atoms.push_back(a.symbol);
  LOGGER(0) << "xyzStructureFile: " << xyz << std::endl;
  LOGGER(0) << "#atoms: " << structure.size() << std::endl;
  for (const auto& a: structure) {
    LOGGER(0) << a.symbol     << ": "
              << a.position.x << ", "
              << a.position.y << ", "
              << a.position.z << std::endl;
  }

  const MovecReader movec(fileName);
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

  /* TODO: reorder
  for (const auto &p: transformation.reorder) { ...  }
  */

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
