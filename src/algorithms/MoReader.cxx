#include <string>
#include <vector>
#include <algorithm>
#include <algorithms/MoReader.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <util/Tensor.hpp>
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
#define IF_GIVEN(_l, ...)                                                      \
  if (isArgumentGiven(_l)) { __VA_ARGS__ }

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(MoReader) {}

int frozenElement(std::string e) {
  if (e == "H" || e == "He") return 0;
  else if (e == "Li" || e == "Be" || e == "B" || e == "C" || e == "N"
           || e == "O" || e == "F" || e == "Ne")
    return 2;
  else if (e == "Na" || e == "Mg" || e == "Al" || e == "Si" || e == "P"
           || e == "S" || e == "Cl" || e == "Ar")
    return 10;
  else throw "I am a quantum chemist. I dont know your crazy element: " + e;
  return -1;
}

std::vector<std::vector<size_t>> findShellIndices(
    const BasisSet &bs,
    const std::vector<std::string> &atoms // atom names in our structure
    ,
    const am::AngularMomentum &am) {
  std::vector<std::vector<size_t>> indices;
  size_t atomBlock(0);

  for (const auto &a : atoms) // e.g. {N , N }
    for (const auto &b : bs) {
      if (b.atom != a) continue;

      // Basis b of atom a in the structure
      for (const auto &shell : b.shells) {
        size_t nbf(am::toInt(shell.am));
        if (shell.am != am) { // e.g. we're looking for D and get S shell
          atomBlock += nbf;
          continue;
        }
        std::vector<size_t> _indices(nbf);
        std::iota(_indices.begin(), _indices.end(), atomBlock);
        indices.push_back(_indices);
        atomBlock += nbf;
      }
    }

  return indices;
}

struct ShellParser {
  const Regex shell_symbol =
                  pars::oneOf({"S", "P", "D", "F", "G", "H", "I", "K"}),
              atom = pars::upper + pars::lower + pars::optional,
              sep = pars::blank + pars::oneOrMore,
              shell = pars::blank
                    + pars::anyOf
                    // number of shells
                    + pars::capture(pars::digit + pars::oneOrMore)
                    // shell sybmol
                    + pars::capture(shell_symbol.s) + pars::blank + pars::anyOf,
              basis_line =
                  pars::capture(atom.s) + sep.s
                  + pars::capture(pars::group(shell.s) + pars::oneOrMore);
  BasisSet parseString(const std::string &s) {
    std::smatch basisMatch, shellMatch;
    std::string S(s);
    BasisSet bs;
    while (std::regex_search(S, basisMatch, basis_line.r)) {
      std::string atomSymbol = basisMatch[1].str(),
                  shellLine = basisMatch[2].str();
      std::cout << basisMatch[0].str() << std::endl;
      std::cout << "shellLine: " << shellLine << std::endl;
      std::vector<Shell> shells;
      while (std::regex_search(shellLine, shellMatch, shell.r)) {
        size_t nShells(std::atoi(shellMatch[1].str().c_str()));
        std::cout << "nshells " << nShells << std::endl;
        am::AngularMomentum am(am::fromString(shellMatch[2].str()));
        for (size_t i(0); i < nShells; i++) {
          shells.push_back({atomSymbol, am, {}});
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
    MoReader::DEFAULT_SCALINGS = {
        {"nwchem",
         {{"DScaling", "1,1,1,-1,1"},
          {"FScaling", "1,1,1,1,-1,1,-1"},
          {"GScaling", "1,1,1,1,1,-1,1,-1,1"},
          {"HScaling", "1,1,1,1,1,1,-1,1,-1,1,-1"},
          {"IScaling", "1,1,1,1,1,1,1,-1,1,-1,1,-1,1"}}},
        {"psi4", {}},
        {"turbomole",
         {{"FScaling", "1, 1, 1, 1, 1, 1,-1"},
          {"GScaling", " 1, 1, 1, 1,-1, 1,-1, 1, 1"},
          {"HScaling", " 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1"}}}};

std::map<std::string, std::map<std::string, std::string>>
    MoReader::DEFAULT_REORDER = {{"nwchem", {}},
                                 {"psi4", {}},
                                 {"turbomole",
                                  {{"DReorder", "3,2,0,1,4"},
                                   {"FReorder", "6,3,2,0,1,4,5"},
                                   {"GReorder", "7,6,3,2,0,1,4,5,8"},
                                   {"HReorder", "10,7,6,3,2,0,1,4,5,8,9"}}}};

std::vector<std::string> MoReader::BACKENDS = {"nwchem", "psi4", "turbomole"};


DEFSPEC(
    MoReader,
    SPEC_IN({"frozenCore", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
            {"backend", SPEC_VALUE("TODO: DOC", std::string)},
            {"basisFile", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
            {"file", SPEC_VALUE("TODO: DOC", std::string)},
            {"shellsFile", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
            {"xyzStructureFile", SPEC_VALUE_DEF("TODO: DOC", std::string, "")}),
    SPEC_OUT(
        {"HoleEigenEnergies", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"OccupationNumbers", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"OrbitalCoefficients", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"ParticleEigenEnergies", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
        {"Spins", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(MoReader) {
  std::vector<std::string> args;
  const std::string fileName(in.get<std::string>("file")),
      xyz(in.get<std::string>("xyzStructureFile", "")),
      basisFile(in.get<std::string>("basisFile", "")),
      shellsFile(in.get<std::string>("shellsFile", "")),
      backend(in.get<std::string>("backend"));
  bool frozenCore(in.get<int64_t>("frozenCore", 0) == 1);

  if (std::find(BACKENDS.begin(), BACKENDS.end(), backend) == BACKENDS.end()) {
    throw "Incorrect backend value: " + backend;
  }

  // make <shell>Scaling
  const auto _s = [&](std::string s) {
    return std::make_pair<am::AngularMomentum, std::vector<double>>(
        am::fromString(s),
        pars::parseVector<double>(this->in.get<std::string>(
            s + "Scaling",
            this->DEFAULT_SCALINGS[backend][s + "Scaling"])));
  };
  // make <shell>Reorder
  const auto _r = [&](std::string s) {
    return std::make_pair<am::AngularMomentum, std::vector<int>>(
        am::fromString(s),
        pars::parseVector<int>(this->in.get<std::string>(
            s + "Reorder",
            this->DEFAULT_REORDER[backend][s + "Reorder"])));
  };
  struct {
    const std::map<am::AngularMomentum, std::vector<double>> scaling;
    const std::map<am::AngularMomentum, std::vector<int>> reorder;
  } transformation = {
      {_s("S"), _s("P"), _s("D"), _s("F"), _s("G"), _s("H"), _s("I"), _s("K")},
      {_r("S"), _r("P"), _r("D"), _r("F"), _r("G"), _r("H"), _r("I"), _r("K")}};

  LOGGER(0) << "NOTE: it only works now for restricted references" << std::endl;
  LOGGER(0) << "file: " << fileName << std::endl;

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
  if (!basis.size()) throw "Shell information not read properly.";
  for (const auto &b : basis) {
    LOGGER(0) << ":: " << b.atom << " :name " << b.name << " "
              << ":#shells" << b.shells.size() << std::endl;
  }

  std::vector<std::string> atoms;
  int frozenElectrons(0);
  if (xyz.size()) {
    const auto structure(pars::XyzParser().parseFile(xyz));
    for (const auto &a : structure) atoms.push_back(a.symbol);
    LOGGER(0) << "xyzStructureFile: " << xyz << std::endl;
    LOGGER(0) << "#atoms: " << structure.size() << std::endl;
    for (const auto &a : structure) {
      LOGGER(0) << a.symbol << ": " << a.position.x << ", " << a.position.y
                << ", " << a.position.z << std::endl;
    }
    if (frozenCore) {
      for (const auto &a : structure) {
        LOGGER(0) << "atoms: " << a.symbol << std::endl;
        frozenElectrons += frozenElement(a.symbol);
      }
      LOGGER(0) << "frozen: " << frozenElectrons << std::endl;
    }
  }

  std::vector<std::vector<double>> mos, eigenvalues, occupations;
  size_t Np;

  if (backend == NWCHEM) {
    const nwchem::MovecReader movec(fileName);
    mos = movec.mos;
    eigenvalues = movec.eigenvalues;
    occupations = movec.occupations;
    Np = movec.Np;
  }
  /*
 * else if (backend == TURBOMOLE) {
   const tmole::MosParser movec(fileName);
   mos = movec.mos;
   eigenvalues = movec.eigenvalues;
   occupations = movec.occupations;
   Np = movec.Np;
 } */

  // For every spin channel, reorder and rescale the mos.
  for (auto &_mos : mos) {

    // do scaling
    for (const auto &p : transformation.scaling) {
      const auto am = am::toInt(p.first);
      const auto &scaling = p.second;
      if (scaling.size()) {
        LOGGER(0) << "doing scaling of am: " << am << std::endl;
        LOGGER(0) << "You should have " << am << " numbers" << std::endl;
        assert(scaling.size() == am);
        auto indicesVector(findShellIndices(basis, atoms, p.first));
        // go through every state (column of mos)
        for (size_t j(0); j < Np; j++) {
          for (auto const &indices : indicesVector) {
            assert(indices.size() == am);
            for (size_t ii(0); ii < am; ii++) {
              assert(indices[ii] < Np);
              _mos[indices[ii] + j * Np] *= scaling[ii];
            } // ii
          }   // indices
        }     // j
        LOGGER(0) << "done with " << am << std::endl;
      }
    }

    // do reorder
    for (const auto &p : transformation.reorder) {
      const auto am = am::toInt(p.first);
      const auto &reorder = p.second;
      if (reorder.size()) {
        LOGGER(0) << "doing reorder of am: " << am << std::endl;
        LOGGER(0) << "You should have " << am << " numbers" << std::endl;
        assert(reorder.size() == am);
        auto indicesVector(findShellIndices(basis, atoms, p.first));
        // go through every state (column of mos)
        for (size_t j(0); j < Np; j++) {
          for (auto const &indices : indicesVector) {
            assert(indices.size() == am);
            // save a backup of the indices thunk to be settable
            std::vector<double> backup(indices.size());
            for (size_t ii(0); ii < am; ii++)
              backup[ii] = _mos[indices[ii] + j * Np];
            for (size_t ii(0); ii < am; ii++) {
              assert(indices[ii] < Np);
              _mos[indices[ii] + j * Np] = backup[reorder[ii]];
            } // ii
          }   // indices
        }     // j
        LOGGER(0) << "done with " << am << std::endl;
      }
    }
  }

  std::vector<int> pp(2, Np), syms(2, NS), o(1), v(1), p(1, Np);
  size_t no(0);
  size_t nv(0);
  std::vector<int64_t> ids;
  const bool unrestricted = mos.size() == 2;
  const int rank_m = int(Sisi4s::world->rank == 0); // rank mask

  std::vector<double> outMos, outEigenvalues, outOccupations, outSpins;
  if (unrestricted) {
    LOGGER(0) << "considering UNRESTRICTED system" << std::endl;
    size_t na(eigenvalues[0].size());
    size_t nb(eigenvalues[1].size());
    size_t npp(na + nb);

    for (size_t b(0); b < occupations[0].size(); b++) {
      if (b < frozenElectrons / 2) continue;
      auto occ(occupations[0][b]);
      if (occ < 0.5) continue;
      no++;
      outOccupations.push_back(occ);
      outSpins.push_back(0.5);
      outEigenvalues.push_back(eigenvalues[0][b]);
      for (size_t jj(0); jj < Np; jj++) outMos.push_back(mos[0][b * Np + jj]);
    }

    for (size_t b(0); b < occupations[1].size(); b++) {
      if (b < frozenElectrons / 2) continue;
      auto occ(occupations[1][b]);
      if (occ < 0.5) continue;
      no++;
      outOccupations.push_back(occ);
      outSpins.push_back(-0.5);
      outEigenvalues.push_back(eigenvalues[1][b]);
      for (size_t jj(0); jj < Np; jj++) outMos.push_back(mos[1][b * Np + jj]);
    }

    for (size_t b(0); b < occupations[0].size(); b++) {
      if (b < frozenElectrons / 2) continue;
      auto occ(occupations[0][b]);
      if (occ > 0.5) continue;
      nv++;
      outOccupations.push_back(occ);
      outSpins.push_back(0.5);
      outEigenvalues.push_back(eigenvalues[0][b]);
      for (size_t jj(0); jj < Np; jj++) outMos.push_back(mos[0][b * Np + jj]);
    }

    for (size_t b(0); b < occupations[1].size(); b++) {
      if (b < frozenElectrons / 2) continue;
      auto occ(occupations[1][b]);
      if (occ > 0.5) continue;
      nv++;
      outOccupations.push_back(occ);
      outSpins.push_back(-0.5);
      outEigenvalues.push_back(eigenvalues[1][b]);
      for (size_t jj(0); jj < Np; jj++) outMos.push_back(mos[1][b * Np + jj]);
    }

    o[0] = no;
    v[0] = nv;
    p[0] = no + nv;
    pp[1] = p[0];
  } else {
    LOGGER(0) << "considering RESTRICTED system" << std::endl;
    for (size_t ii(0); ii < occupations[0].size(); ii++) {
      if (ii < frozenElectrons / 2) continue;
      auto occ(occupations[0][ii]);
      if (occ > 0.5) {
        no++;
      } else {
        nv++;
      }
      outEigenvalues.push_back(eigenvalues[0][ii]);
      outOccupations.push_back(occ);
      for (size_t jj(0); jj < Np; jj++) outMos.push_back(mos[0][ii * Np + jj]);
    }
    o[0] = no;
    v[0] = nv;
    p[0] = no + nv;
    pp[1] = p[0];
  }

  IF_GIVEN(
      "HoleEigenEnergies", ids.resize(rank_m * o[0]);
      std::iota(ids.begin(), ids.end(), 0);
      auto epsi(new Tensor<double>(1, o.data(), syms.data(), *Sisi4s::world));
      epsi->write(ids.size(), ids.data(), outEigenvalues.data());
      out.set<Tensor<double> *>("HoleEigenEnergies", epsi);)

  IF_GIVEN(
      "ParticleEigenEnergies", ids.resize(rank_m * v[0]);
      std::iota(ids.begin(), ids.end(), 0);
      auto epsa(new Tensor<double>(1, v.data(), syms.data(), *Sisi4s::world));
      epsa->write(ids.size(), ids.data(), outEigenvalues.data() + o[0]);
      out.set<Tensor<double> *>("ParticleEigenEnergies", epsa);)

  IF_GIVEN(
      "OccupationNumbers", ids.resize(rank_m * p[0]);
      std::iota(ids.begin(), ids.end(), 0);
      auto os(new Tensor<double>(1, p.data(), syms.data(), *Sisi4s::world));
      os->write(ids.size(), ids.data(), outOccupations.data());
      out.set<Tensor<double> *>("OccupationNumbers", os);)
  if (unrestricted) {
    ids.resize(rank_m * p[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto spin(new Tensor<double>(1, p.data(), syms.data(), *Sisi4s::world));
    spin->write(ids.size(), ids.data(), outSpins.data());
    out.set<Tensor<double> *>("Spins", spin);
  }

  if (isArgumentGiven("OrbitalCoefficients")) {
    ids.resize(rank_m * pp[0] * pp[1]);
    std::iota(ids.begin(), ids.end(), 0);
    auto coef(new Tensor<double>(2, pp.data(), syms.data(), *Sisi4s::world));
    coef->write(ids.size(), ids.data(), outMos.data());
    out.set<Tensor<double> *>("OrbitalCoefficients", coef);
  }
}

SISI_DOC(*Documentation

             This module contains convenience routines in order to read
                 molecular orbitals from the following codes
         :

         -NWCHEM - PSI4
             - TURBOMOLE

                   * END)
