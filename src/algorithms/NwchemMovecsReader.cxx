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

void NwchemMovecsReader::run() {
  const std::string fileName(getTextArgument("file"));
  const int No(getIntegerArgument("No"));

  LOGGER(0) << "NOTE: it only works now for restricted references" << std::endl;
  LOGGER(0) << "file: " << fileName << std::endl;
  LOGGER(0) << "No: " << No << std::endl;

  std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
  size_t nsets, Np;
  std::vector<size_t> nmo;
  std::vector<double> occupations, eigenvalues, mos;
  size_t iset = 0; // This allows to skip the first set, dont need this yet

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
      for (size_t j(0); j < mos_buff.size(); j++) mos[j + Np*i] = mos_buff[j];
    }

  }

  file.close();

  std::vector<int> pp(2, Np), syms(2, NS), o(1, No), v(1, Np - No);
  std::vector<int64_t> ids;


  IF_GIVEN("HoleEigenEnergies",
    ids.resize(o[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsi(new CTF::Tensor<double>(1, o.data(), syms.data(), *Cc4s::world));
    epsi->write(ids.size(), ids.data(), eigenvalues.data());
    allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  )

  IF_GIVEN("ParticleEigenEnergies",
    ids.resize(v[0]);
    std::iota(ids.begin(), ids.end(), 0);
    auto epsa(new CTF::Tensor<double>(1, v.data(), syms.data(), *Cc4s::world));
    epsa->write(ids.size(), ids.data(), eigenvalues.data() + No);
    allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);
  )

  IF_GIVEN("OccupationNumbers",
    ids.resize(Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto os(new CTF::Tensor<double>(1, pp.data(), syms.data(), *Cc4s::world));
    os->write(ids.size(), ids.data(), occupations.data());
    allocatedTensorArgument<double>("OccupationNumbers", os);
  )

  IF_GIVEN("OrbitalCoefficients",
    ids.resize(Np*Np);
    std::iota(ids.begin(), ids.end(), 0);
    auto coef(new CTF::Tensor<double>(2, pp.data(), syms.data(), *Cc4s::world));
    coef->write(ids.size(), ids.data(), mos.data());
    allocatedTensorArgument<double>("OrbitalCoefficients", coef);
  )

}
