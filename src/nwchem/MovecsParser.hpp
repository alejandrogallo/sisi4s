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
#include <boost/lexical_cast.hpp>

#define _L(_l) LOG(_l, "MovecReader")

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
using namespace pars;
namespace nwchem {

struct MovecReader {
  size_t nsets, Np;
  std::vector<size_t> nmo;
  std::vector<double> occupations, eigenvalues, mos;
  size_t iset = 0; // This allows to skip the first set, dont need this yet
  MovecReader(const std::string &fileName) {
    std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    if (!file.good())    throw "File " + fileName + " not found";
    if (!file.is_open()) throw "File IO error: " + fileName;
    readFortranChunk<char>(file); // convergence info
    readFortranChunk<char>(file); // scftype
    readFortranChunk<char>(file); // lentit
    { const auto _rfile(readFortranChunk<char>(file).data());
      _L(0) << "title: " << _rfile << std::endl;
    }
    readFortranChunk<char>(file).data(); // lenbas
    { const auto _rbasis(readFortranChunk<char>(file).data());
      _L(0) << "basis: " << _rbasis << std::endl;
    }

    nsets = *readFortranChunk<int64_t>(file).data(); // nsets
    Np = *readFortranChunk<int64_t>(file).data();   // Np

    _L(0) << "nsets: " << nsets << std::endl;
    _L(0) << "Np: " << Np << std::endl;

    // read nmos
    for (size_t i=0; i<nsets; i++) {
      int64_t buf = *readFortranChunk<int64_t>(file).data();
      nmo.push_back(buf);
      _L(0) << "nmo[" << i << "]: " << nmo[0] << std::endl;
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
      _L(0) << "#occupations read: " << occupations.size() << std::endl;
      assert(occupations.size() == Np);

      // Eigenvaluess [read(unitno) (evals(j), j=1,Np)]
      eigenvalues = readFortranChunk<double>(file);
      _L(0) << "#eigenvalues read: " <<  eigenvalues.size() << std::endl;
      assert(eigenvalues.size() == Np);

      mos.resize(Np * Np);
      for (size_t i(0); i < nmo[iset]; i++) { // do i = 1, nmo(iset)
        const std::vector<double> mos_buff(readFortranChunk<double>(file));
        assert(mos_buff.size() == Np);
        for (size_t j(0); j < mos_buff.size(); j++)
          mos[j + Np*i] = mos_buff[j];
      }
    }

    file.close();

  }
};


struct TmoleMosReader {
  size_t Np;
  std::vector<size_t> Nsaos;
  std::vector<double> occupations, eigenvalues, mos;
  std::smatch match;
  const pars::Regex sep       = blank + oneOrMore,
                    number    = digit + oneOrMore,
                    twodigits = digit + digit,
                    expnumber = anyChar + "." + digit + oneOrMore +
                                "D" + anyChar + twodigits.s,
                    symm      = alnum + oneOrMore,
                    comment   = bof + "#" + print + anyOf,
                    unper     = bof + anyChar + "scfmo" + print + anyOf,
                    eigenval  = bof + sep.s + capture(number.s) + sep.s +
                                capture(symm.s) + sep.s + "eigenvalue=" +
                                capture(expnumber.s) + sep.s + "nsaos=" +
                                capture(number.s) + print + anyOf,
                    mosline1   = capture(expnumber.s),
                    mosline2   = capture(expnumber.s) + capture(expnumber.s),
                    mosline3   = capture(expnumber.s) + capture(expnumber.s)
                               + capture(expnumber.s),
                    mosline4   = capture(expnumber.s) + capture(expnumber.s)
                               + capture(expnumber.s) + capture(expnumber.s);



  size_t iset = 0; // This allows to skip the first set, dont need this yet
  TmoleMosReader(const std::string &fileName) {
    std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    std::vector<double> mos_buf;
    std::string s, e;
    size_t buf_size;
    if (!file.good())    throw "File " + fileName + " not found";
    if (!file.is_open()) throw "File IO error: " + fileName;
    std::string line;
    while (std::getline(file, line)) {
      if ( std::regex_match(line, match, comment.r)
        || std::regex_match(line, match, unper.r)  ) continue;
      if ( std::regex_match(line, match, eigenval.r)) {
        e = std::string(match[3]);
        mos_buf.resize(0);
        buf_size = atoi(std::string(match[4]).c_str());
        Nsaos.push_back(buf_size);
        std::cout << buf_size << std::endl;
      }

      if ( std::regex_match(line, match, mosline1.r))
        for (int i(1); i < match.size(); i++) {
          s = std::string(match[i]);
          std::replace( s.begin(), s.end(), 'D', 'e');
          mos_buf.push_back(boost::lexical_cast<double>(s.c_str()));
        }
      if ( std::regex_match(line, match, mosline2.r))
			  for (int i(1); i < match.size(); i++) {
          s = std::string(match[i]);
          std::replace( s.begin(), s.end(), 'D', 'e');
          mos_buf.push_back(boost::lexical_cast<double>(s.c_str()));
        }
      if ( std::regex_match(line, match, mosline3.r))
        for (int i(1); i < match.size(); i++) {
          s = std::string(match[i]);
          std::replace( s.begin(), s.end(), 'D', 'e');
          mos_buf.push_back(boost::lexical_cast<double>(s.c_str()));
        }
      if ( std::regex_match(line, match, mosline4.r))
        for (int i(1); i < match.size(); i++) {
          s = std::string(match[i]);
          std::replace( s.begin(), s.end(), 'D', 'e');
          mos_buf.push_back(boost::lexical_cast<double>(s.c_str()));
        }
      if (mos_buf.size() > 0 && mos_buf.size() == buf_size) {
        std::replace(e.begin(), e.end(), 'D', 'e');
        eigenvalues.push_back(boost::lexical_cast<double>(e.c_str()));
        for (int64_t i(0); i < mos_buf.size(); i++){
           mos.push_back(mos_buf[i]);
        }
        mos_buf.resize(0);
      }
    }
    Np = eigenvalues.size();

    file.close();

  }
};


}

#endif
