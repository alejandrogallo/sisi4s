#include <algorithms/FcidumpReader.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <ctf.hpp>
#include <regex>
#include <algorithm>
#include <numeric>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(FcidumpReader);

FcidumpReader::FcidumpHeader
FcidumpReader::parseHeader(const std::string &filePath) {
  FcidumpReader::FcidumpHeader header{.norb=0, .nelec=0, .uhf=0, .ms2=0};
  std::regex
    rnorb{"NORB\\s*=\\s*([0-9]+)"},
    rnelec{"NELEC\\s*=\\s*([0-9]+)"},
    ruhf{"UHF\\s*=\\s*([01])"},
    rms2{"MS2\\s*=\\s*([0-9]+)"},
    rend{"^\\s*([/]|&END|\\$)\\s*$"};
  std::string line;
  std::ifstream file(filePath);
  if (file.is_open()) {
    while(std::getline(file, line)) {
      std::smatch matches;
      std::regex_search(line, matches, rnorb);  // parse norb
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

struct IntegralParser {
  static const size_t columns{4};
  std::string name;
  std::regex regex;
  std::vector<double> values;
  std::vector<size_t> indices;
  std::vector<int> lens;
  std::vector<int> syms;
  int No, Nv;
  size_t dimension;
  IntegralParser(std::string name_, const FcidumpReader::FcidumpHeader &header)
      :name(name_) {
    if (!std::regex_match(name, std::regex{"^[hp]*$"})) {
      throw new EXCEPTION("Name should be a combination of [hp] or empty");
    }
    syms.insert(syms.begin(), name.size(), NS);
    No = (header.uhf == 1 ? 1 : 0.5) * header.nelec;
    Nv = (header.uhf == 1 ? 2 : 1) * header.norb - No;
    std::for_each(name.begin(), name.end(), // build lens
      [&](char a){ lens.push_back(a == 'h' ? No : Nv); }
    );
    dimension = std::accumulate(
      lens.begin(), lens.end(), 1, std::multiplies<int>());
    std::string regex_str{"^\\s*(\\S+)"};
    std::for_each(
      lens.begin(), lens.end(), [&](int){ regex_str += "\\s+([1-9][0-9]*)"; });
    for (size_t i(0); i<(columns - lens.size()); ++i) {
      regex_str += "\\s+0";
    }
    regex_str += "\\s*$";
    regex = std::regex{regex_str};
    //LOG(0, "FcidumpReader") << regex_str << std::endl;
  }
  bool match(const std::string &line) {
    std::smatch matches;
    std::regex_match(line, matches, regex);
    if (matches.size()) {
      std::vector<int> gIndices(lens.size());
      // matches = {line, value, index1, index2, ...}
      for(unsigned int i(2); i<matches.size(); i++) {
        int k =
          std::atoi(std::string{matches[i]}.c_str())
        ;
        if ((k <= No && name[i-2] == 'p') || (k > No && name[i-2] == 'h')){
          return false;
        }
        gIndices[i-2] = k;
      }
      LOG(1, "FcidumpReader") << name << ": " << line << std::endl;
      std::vector<int> gIndexLens(lens.size());
      {// build up gIndexLens (1, N_1, N_1 * N_2, ..., N_1 *...* N_n-1)
        std::vector<int> tempLens(lens);
        tempLens.pop_back(); tempLens.insert(tempLens.begin(), 1);
        std::partial_sum(
          tempLens.begin(), tempLens.end(),
          gIndexLens.begin(), std::multiplies<int>());
      }
      values.push_back(std::atof(std::string{matches[1]}.c_str()));
      indices.push_back(
        std::inner_product(
          gIndexLens.begin(), gIndexLens.end(), gIndices.begin(), 0));
      return true;
    } else {
      return false;
    }
  }
};

void FcidumpReader::run() {
  auto filePath(getTextArgument("file", "FCIDUMP"));
  LOG(0, "FcidumpReader") << "Reading fcidump from " << filePath << std::endl;
  FcidumpReader::FcidumpHeader header(parseHeader(filePath));
  LOG(0, "FcidumpReader") << "NELEC = " << header.nelec << std::endl;
  LOG(0, "FcidumpReader") << "NORB  = " << header.norb << std::endl;
  LOG(0, "FcidumpReader") << "MS2   = " << header.ms2 << std::endl;
  LOG(0, "FcidumpReader") << "UHF   = " << header.uhf << std::endl;
  int No((header.uhf == 1 ? 1 : 0.5) * header.nelec);
  int Nv((header.uhf == 1 ? 2 : 1) * header.norb - No);
  LOG(0, "FcidumpReader") << "No    = " << No << std::endl;
  LOG(0, "FcidumpReader") << "Nv    = " << Nv << std::endl;

  const std::vector<std::string> integralNames{
    "hh", "pp", "hp", "ph",
    "hhhh", "hhhp", "hhph", "hhpp", "hphh", "hphp", "hpph", "hppp",
    "phhh", "phhp", "phph", "phpp", "pphh", "pphp", "ppph", "pppp" };
  std::vector<IntegralParser> integralParsers;

  for (const auto& name: integralNames) {
    if (isArgumentGiven(name)) {
      integralParsers.push_back(IntegralParser(name, header));
      LOG(0, "FcidumpReader") << "Parsing " << name << std::endl;
    }
  }

  std::ifstream file(filePath);
  std::string line;
  if (file.is_open()) {
    while(std::getline(file, line)) {
      for (auto it=integralParsers.begin(); it != integralParsers.end(); ++it){
        if ((*it).match(line)){
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
}
