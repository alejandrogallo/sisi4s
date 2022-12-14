#ifndef _MOVECS_PARSER_HEADE
#define _MOVECS_PARSER_HEADE

#include <regex>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <assert.h>
#include <util/Log.hpp>
#include <util/Parsing.hpp>

#define _L(_l) LOG(_l, "MovecReader")

template <typename F>
std::vector<F> readFortranChunk(std::fstream &file) {
  uint32_t fortran_section, copy;
  std::vector<F> buffer;

  file.read((char *)&fortran_section, sizeof(int32_t));
  copy = fortran_section;
  buffer.resize(fortran_section / sizeof(F));
  file.read((char *)buffer.data(), fortran_section);
  file.read((char *)&fortran_section, sizeof(int32_t));
  assert(copy == fortran_section);

  return buffer;
}
using namespace pars;
namespace nwchem {

struct MovecReader {
  size_t nsets, Np;
  std::vector<int64_t> nmo;
  std::vector<std::vector<double>> occupations, eigenvalues, mos;
  size_t iset = 0; // This allows to skip the first set, dont need this yet
  MovecReader(const std::string &fileName) {
    std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    if (!file.good()) throw "File " + fileName + " not found";
    if (!file.is_open()) throw "File IO error: " + fileName;
    readFortranChunk<char>(file); // convergence info
    readFortranChunk<char>(file); // scftype
    readFortranChunk<char>(file); // lentit
    {
      const auto _rfile(readFortranChunk<char>(file).data());
      _L(0) << "title: " << _rfile << std::endl;
    }
    readFortranChunk<char>(file).data(); // lenbas
    {
      const auto _rbasis(readFortranChunk<char>(file).data());
      _L(0) << "basis: " << _rbasis << std::endl;
    }

    nsets = *readFortranChunk<int64_t>(file).data(); // nsets
    Np = *readFortranChunk<int64_t>(file).data();    // Np

    _L(0) << "nsets: " << nsets << std::endl;
    _L(0) << "Np: " << Np << std::endl;

    // resize occupations and such
    occupations.resize(nsets);
    eigenvalues.resize(nsets);
    mos.resize(nsets);

    // read nmos
    nmo = readFortranChunk<int64_t>(file);
    for (size_t i = 0; i < nsets; i++) {
      _L(0) << "nmo[" << i << "]: " << nmo[i] << std::endl;
    }

    // loop which reject MOs from the set
    for (size_t i = 0; i < iset; i++) {
      readFortranChunk<char>(file);
      readFortranChunk<char>(file);
      for (uint64_t j = 0; j < nmo[i]; j++) { // do i = 1, nmo(jset)
                                              // read(unitno)
        readFortranChunk<char>(file);
      }
    }

    for (size_t s(0); s < nsets; s++) { // do s = 1, nsets

      // Read occupation numbers
      occupations[s] = readFortranChunk<double>(file);
      _L(0) << "#occupations read: " << occupations[s].size() << std::endl;
      assert(occupations[s].size() == Np);
      // Eigenvaluess [read(unitno) (evals(j), j=1,Np)]
      eigenvalues[s] = readFortranChunk<double>(file);
      _L(0) << "#eigenvalues read: " << eigenvalues[s].size() << std::endl;
      assert(eigenvalues[s].size() == Np);
      mos[s].resize(Np * Np);
      for (size_t i(0); i < nmo[iset]; i++) { // do i = 1, nmo(iset)
        const std::vector<double> mos_buff(readFortranChunk<double>(file));
        assert(mos_buff.size() == Np);
        for (size_t j(0); j < mos_buff.size(); j++)
          mos[s][j + Np * i] = mos_buff[j];
      }
    }

    file.close();
  }
};

} // namespace nwchem

#endif
