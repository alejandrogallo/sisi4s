#include <algorithms/FcidumpReader.hpp>
#include <util/Tensor.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Exception.hpp>
#include <util/Integrals.hpp>
#include <fstream>
#include <regex>
#include <algorithm>
#include <numeric>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(FcidumpReader) {}

FcidumpReader::FcidumpHeader
FcidumpReader::parseHeader(const std::string &filePath) {
  FcidumpReader::FcidumpHeader header{.norb = 0,
                                      .nelec = 0,
                                      .uhf = 0,
                                      .ms2 = 0};
  std::regex rnorb{"NORB\\s*=\\s*([0-9]+)"}, rnelec{"NELEC\\s*=\\s*([0-9]+)"},
      ruhf{"UHF\\s*=\\s*([01])"}, rms2{"MS2\\s*=\\s*([0-9]+)"},
      rend{"^\\s*([/]|&END|\\$)\\s*$"};
  std::string line;
  std::ifstream file(filePath);
  if (file.is_open()) {
    while (std::getline(file, line)) {
      std::smatch matches;
      std::regex_search(line, matches, rnorb); // parse norb
      if (matches.size()) header.norb = atoi(std::string{matches[1]}.c_str());
      std::regex_search(line, matches, rnelec); // parse nelec
      if (matches.size()) header.nelec = atoi(std::string{matches[1]}.c_str());
      std::regex_search(line, matches, ruhf); // parse uhf
      if (matches.size()) header.uhf = atoi(std::string{matches[1]}.c_str());
      std::regex_search(line, matches, rms2); // parse ms2
      if (matches.size()) header.ms2 = atoi(std::string{matches[1]}.c_str());
      std::regex_match(line, matches, rend); // parse end of header
      if (matches.size()) {
        LOG(0, "FcidumpReader") << "Breaking at line " << line << std::endl;
        break;
      }
      LOG(0, "FcidumpReader") << line << std::endl;
    }
  }
  file.close();
  return header;
}

int indexToLenght(const char a, const int No, const int Nv) {
  if (a == 'h') {
    return No;
  } else if (a == 'p') {
    return Nv;
  } else {
    return No + Nv;
  }
}

struct IntegralParser {
  // how many index columns are there in the fcidump
  static const size_t index_columns{4};
  // e.g. hh
  std::string name;
  // e.g. hh in chemist notation
  std::string chemistName;
  // regex for a single line
  std::regex regex;
  // actual values of the tensor
  std::vector<double> values;
  // actual indices
  std::vector<int64_t> indices;
  // ctf lens
  std::vector<int> lens;
  // ctf syms
  std::vector<int> syms;
  bool uhf;
  int No, Nv;
  // general size of the tensor
  size_t dimension;
  IntegralParser(std::string name_, const FcidumpReader::FcidumpHeader &header)
      : name(name_) {

    if (!std::regex_match(name, std::regex{"^[hpt]*$"})) {
      throw new EXCEPTION("Name should be a combination of [hpt] or empty");
    }

    // create the chemistName of this integral
    chemistName = (name.size() == 4) ? permuteIndices(name, 1, 2) : name;

    uhf = header.uhf;
    No = (uhf == 1 ? 1 : 0.5) * header.nelec;
    Nv = (uhf == 1 ? 2 : 1) * header.norb - No;
    // build syms
    syms.insert(syms.begin(), name.size(), NS);
    // build lens
    std::for_each(name.begin(), name.end(), [&](const char &a) {
      lens.push_back(indexToLenght(a, No, Nv));
    });
    // build dimensio
    dimension =
        std::accumulate(lens.begin(), lens.end(), 1, std::multiplies<int>());
    // start to build regex for the numbers
    // it starts with the numerical value of the integral
    std::string regex_str{"^\\s*(\\S+)"};
    // and then so many indices as lens we have
    std::for_each(lens.begin(), lens.end(), [&](int) {
      regex_str += "\\s+([1-9][0-9]*)";
    });
    // if lens is less that the number of columns available, the rest should
    // be zeros
    for (size_t i(0); i < (index_columns - lens.size()); ++i) {
      regex_str += "\\s+0";
    }
    // end with possible padding zeros
    regex_str += "\\s*$";
    regex = std::regex{regex_str};
  }

  bool match(const std::string &line) {
    std::smatch matches;
    std::regex_match(line, matches, regex);
    std::vector<int> gIndices(lens.size());
    std::vector<int> gIndexLens(lens.size());
    const int spin_uhf(uhf ? 1 : 2);

    // matches = {line, value, index1, index2, ...}
    if (matches.size() == 0) return false;

    // check if the line relates to this integral
    for (unsigned int i(2); i < matches.size(); i++) {
      int k = std::atoi(std::string{matches[i]}.c_str());
      // what is the corresponding index in our integrals, H or P?
      const char _HorPorT(chemistName[i - 2]);
      // if the index is not what we're expecting then return false
      if ((k <= No && _HorPorT == 'p') || (k > No && _HorPorT == 'h'))
        return false;
      if (_HorPorT == 'p') {
        gIndices[i - 2] = k - No - 1;
      } else if (_HorPorT == 'h') {
        gIndices[i - 2] = k - 1;
      } else {
        // just get the pure index if
        gIndices[i - 2] = k - 1;
      }
    }
    // LOG(1, "FcidumpReader") << name << ":(" << chemistName << "):"
    //<< line << std::endl;

    // this will be in physics notation
    { // build up gIndexLens (1, N_1, N_1 * N_2, ..., N_1 *...* N_n-1)
      // copy lens into tempLens
      std::vector<int> tempLens(lens);
      // remove the last lens, the last element
      tempLens.pop_back();
      // insert at the beginning a 1
      tempLens.insert(tempLens.begin(), 1);
      // create the list with the partial products of the elements
      // in tempLens, and store them in gIndexLens.begin()
      std::partial_sum(tempLens.begin(),
                       tempLens.end(),
                       gIndexLens.begin(),
                       std::multiplies<int>());
    }
    // gIndices was read in chemist notation since it comes from a chemistName
    // so we have to change it back to physics notation to store the index
    // correctly
    if (gIndices.size() == 4) gIndices = permuteIndices(gIndices, 1, 2);
    values.push_back(std::atof(std::string{matches[1]}.c_str()));
    indices.push_back(std::inner_product(gIndexLens.begin(),
                                         gIndexLens.end(),
                                         gIndices.begin(),
                                         0));
    return true;
  }

  Tensor<double> *allocateTensor() {
    const int rank_m = int(Sisi4s::world->rank == 0); // rank mask
    auto t(new Tensor<double>(lens.size(),
                              lens.data(),
                              syms.data(),
                              *Sisi4s::world));
    t->write(rank_m * indices.size(), indices.data(), values.data());
    return t;
  }
};


DEFSPEC(FcidumpReader,
        SPEC_IN({"nelec", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},
                {"file", SPEC_VALUE_DEF("TODO: DOC", std::string, "FCIDUMP")}),
        SPEC_OUT({name, SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(FcidumpReader) {
  const auto filePath(in.get<std::string>("file", "FCIDUMP"));
  // override the header of the fcidump
  const int nelec(in.get<int64_t>("nelec", -1));
  FcidumpReader::FcidumpHeader header(parseHeader(filePath));
  if (nelec != -1) header.nelec = nelec;
  const int No((header.uhf == 1 ? 1 : 0.5) * header.nelec);
  const int Nv((header.uhf == 1 ? 2 : 1) * header.norb - No);

  if (header.uhf) { throw new EXCEPTION("UHF not implemented"); }

  LOG(0, "FcidumpReader") << "Fcidump = " << filePath << std::endl;
  LOG(0, "FcidumpReader") << "NELEC   = " << header.nelec << std::endl;
  LOG(0, "FcidumpReader") << "NORB    = " << header.norb << std::endl;
  LOG(0, "FcidumpReader") << "MS2     = " << header.ms2 << std::endl;
  LOG(0, "FcidumpReader") << "UHF     = " << header.uhf << std::endl;
  LOG(0, "FcidumpReader") << "No      = " << No << std::endl;
  LOG(0, "FcidumpReader") << "Nv      = " << Nv << std::endl;

  const std::vector<IntegralInfo> twoBody({
      {"tttt", {NP, NP, NP, NP}, "pqrs"},
      {"hhhh", {NO, NO, NO, NO}, "ijkl"},
      {"hhhp", {NO, NO, NO, NV}, "ijka"},
      {"hhph", {NO, NO, NV, NO}, "ijak"},
      {"hhpp", {NO, NO, NV, NV}, "ijab"},
      {"hphh", {NO, NV, NO, NO}, "iajk"},
      {"hphp", {NO, NV, NO, NV}, "iajb"},
      {"hpph", {NO, NV, NV, NO}, "iabj"},
      {"hppp", {NO, NV, NV, NV}, "iabc"},
      {"phhh", {NV, NO, NO, NO}, "aijk"},
      {"phhp", {NV, NO, NO, NV}, "aijb"},
      {"phph", {NV, NO, NV, NO}, "aibj"},
      {"phpp", {NV, NO, NV, NV}, "aibc"},
      {"pphh", {NV, NV, NO, NO}, "abij"},
      {"pphp", {NV, NV, NO, NV}, "abic"},
      {"ppph", {NV, NV, NV, NO}, "abci"},
      {"pppp", {NV, NV, NV, NV}, "abcd"},
  });

  const std::vector<std::string> integralNames{
      "tt",   "tttt", "hh",   "pp",   "hp",   "ph",   "hhhh", "hhhp",
      "hhph", "hhpp", "hphh", "hphp", "hpph", "hppp", "phhh", "phhp",
      "phph", "phpp", "pphh", "pphp", "ppph", "pppp"};
  std::vector<IntegralParser> integralParsers;

  for (const auto &name : integralNames) {
    if (isArgumentGiven(name)) {
      integralParsers.push_back(IntegralParser(name, header));
      LOG(0, "FcidumpReader") << "Parsing " << name << std::endl;
    }
  }

  std::ifstream file(filePath);
  std::string line;
  if (file.is_open()) {
    while (std::getline(file, line)) {
      for (auto it = integralParsers.begin(); it != integralParsers.end();
           ++it) {
        if ((*it).match(line)) {
          // Swap position of the first element in the parsers because
          // it's highly probable that the line in the fcidump also referes
          // to this integral, for large files and many integrals this will
          // be beneficial
          std::iter_swap(integralParsers.begin(), it);
          break;
        }
      }
    }
  }
  file.close();
  LOG(0, "FcidumpReader") << "parsing done" << std::endl;

  for (auto &parser : integralParsers) {
    const auto &name(parser.name);
    if (isArgumentGiven(name)) {
      LOG(0, "FcidumpReader") << "Exporting: " << parser.name << std::endl;
      out.set<Tensor<double> *>(name, parser.allocateTensor());
    }
  }
}
