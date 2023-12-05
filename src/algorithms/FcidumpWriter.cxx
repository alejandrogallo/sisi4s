#include <algorithms/FcidumpReader.hpp>
#include <algorithms/FcidumpWriter.hpp>
#include <util/Tensor.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Exception.hpp>
#include <util/Integrals.hpp>
#include <fstream>
#include <regex>
#include <algorithm>
#include <numeric>
#include <ostream>
#include <cstdio>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(FcidumpWriter);

IMPLEMENT_EMPTY_DRYRUN(FcidumpWriter) {}

constexpr int NxUndefined = -1;
struct TensorInfo {
  const std::string name;
  const std::vector<Index> indices;
};

inline std::ostream &operator<<(std::ostream &s,
                                const FcidumpReader::FcidumpHeader &h) {
  return s << "&FCI" << std::endl
           << " NORB=" << h.norb << std::endl
           << " NELEC=" << h.nelec << std::endl
           << " MS2=" << h.ms2 << std::endl
           << "/" << std::endl;
}

void FcidumpWriter::run() {
  const auto filePath(getTextArgument("file", "FCIDUMP"));
  size_t uhf(getIntegerArgument("uhf", 0));
  size_t ms2(getIntegerArgument("ms2", 0));
  const double threshold(getRealArgument("threshold", 1e-6));
  int No(NxUndefined);
  int Nv(NxUndefined);
  FcidumpReader::FcidumpHeader header;

  const std::vector<TensorInfo> allIntegrals({
      {"PPPP", {NP, NP, NP, NP}}, {"hhhh", {NO, NO, NO, NO}},
      {"hhhp", {NO, NO, NO, NV}}, {"hhph", {NO, NO, NV, NO}},
      {"hhpp", {NO, NO, NV, NV}}, {"hphh", {NO, NV, NO, NO}},
      {"hphp", {NO, NV, NO, NV}}, {"hpph", {NO, NV, NV, NO}},
      {"hppp", {NO, NV, NV, NV}}, {"phhh", {NV, NO, NO, NO}},
      {"phhp", {NV, NO, NO, NV}}, {"phph", {NV, NO, NV, NO}},
      {"phpp", {NV, NO, NV, NV}}, {"pphh", {NV, NV, NO, NO}},
      {"pphp", {NV, NV, NO, NV}}, {"ppph", {NV, NV, NV, NO}},
      {"pppp", {NV, NV, NV, NV}}, {"hh", {NO, NO}},
      {"pp", {NV, NV}},           {"hp", {NO, NV}},
      {"ph", {NV, NO}},
  });

  // find out No and Nv
  for (const auto &integral : allIntegrals) {
    auto &is(integral.indices);
    if (isArgumentGiven(integral.name)) {
      auto it(std::find(is.begin(), is.end(), NO));
      if (No == NxUndefined && it != is.end()) {
        auto tensor(getTensorArgument<double>(integral.name));
        No = tensor->lens[it - is.begin()];
        LOG(0, "FcidumpWriter") << _FORMAT("Detected No: %d from integral %s\n",
                                           No,
                                           integral.name.c_str());
      }
      it = std::find(is.begin(), is.end(), NV);
      if (Nv == NxUndefined && it != is.end()) {
        auto tensor(getTensorArgument<double>(integral.name));
        Nv = tensor->lens[it - is.begin()];
        LOG(0, "FcidumpWriter") << _FORMAT("Detected Nv: %d from integral %s\n",
                                           Nv,
                                           integral.name.c_str());
      }
    }
  }

  // If I don't know No or Nv, panic
  if (No == NxUndefined || Nv == NxUndefined) {
    throw new EXCEPTION("No or Nv could not be parsed!");
  }

  header.uhf = uhf;
  header.nelec = uhf == 1 ? No : No * 2;
  header.norb = No + Nv;
  header.ms2 = ms2;

  LOG(0, "FcidumpWriter") << "Fcidump   = " << filePath << std::endl;
  LOG(0, "FcidumpWriter") << "threshold = " << threshold << std::endl;
  LOG(0, "FcidumpWriter") << "NELEC     = " << header.nelec << std::endl;
  LOG(0, "FcidumpWriter") << "NORB      = " << header.norb << std::endl;
  LOG(0, "FcidumpWriter") << "MS2       = " << header.ms2 << std::endl;
  LOG(0, "FcidumpWriter") << "UHF       = " << header.uhf << std::endl;
  LOG(0, "FcidumpWriter") << "No        = " << No << std::endl;
  LOG(0, "FcidumpWriter") << "Nv        = " << Nv << std::endl;

  // TODO: do this really well, now it's kind of a hack
  // TODO: also write the header into the file
  for (const auto &integral : allIntegrals) {
    if (!isArgumentGiven(integral.name)) continue;
    auto tensor(getTensorArgument<double>(integral.name));
    FILE *fd;
    fd = std::fopen(_FORMAT("%s.fcidump", integral.name.c_str()).c_str(), "w+");
    tensor->print(fd);
    std::fclose(fd);
  }

  std::ifstream file(filePath);
  std::string line;
  if (file.is_open()) {}
  file.close();
}
