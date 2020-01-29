#include <algorithms/FcidumpReader.hpp>
#include <ctf.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <util/Exception.hpp>
#include <util/Integrals.hpp>
#include <fstream>
#include <regex>
#include <algorithm>
#include <numeric>

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
  static const size_t index_columns{4};
  std::string name;
  std::string chemistName;
  std::regex regex;
  std::vector<double> values;
  std::vector<int64_t> indices;
  std::vector<int> lens;
  std::vector<int> syms;
  bool uhf;
  int No, Nv;
  size_t dimension;
  IntegralParser(std::string name_, const FcidumpReader::FcidumpHeader &header)
      :name(name_) {

    if (!std::regex_match(name, std::regex{"^[hp]*$"})) {
      throw new EXCEPTION("Name should be a combination of [hp] or empty");
    }

    // create the chemistName of this integral
    chemistName = (name.size() == 4) ?  permuteIndices(name, 1, 2) : name;

    uhf = header.uhf;
    No = (uhf == 1 ? 1 : 0.5) * header.nelec;
    Nv = (uhf == 1 ? 2 : 1) * header.norb - No;
    // build syms
    syms.insert(syms.begin(), name.size(), NS);
    // build lens
    std::for_each(name.begin(), name.end(),
                  [&](const char& a){ lens.push_back(a == 'h' ? No : Nv); });
    // build dimensio
    dimension = std::accumulate(lens.begin(), lens.end(), 1,
                                std::multiplies<int>());
    // start to build regex for the numbers
    // it starts with the numerical value of the integral
    std::string regex_str{"^\\s*(\\S+)"};
    // and then so many indices as lens we have
    std::for_each(lens.begin(), lens.end(),
                  [&](int){ regex_str += "\\s+([1-9][0-9]*)"; });
    // if lens is less that the number of columns available, the rest should
    // be zeros
    for (size_t i(0); i<(index_columns - lens.size()); ++i) {
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
    for(unsigned int i(2); i<matches.size(); i++) {
      int k = std::atoi(std::string{matches[i]}.c_str());
      // what is the corresponding index in our integrals, H or P?
      const char _HorP(chemistName[i-2]);
      // if the index is not what we're expecting then return false
      if ((k <= No && _HorP == 'p') || (k > No && _HorP == 'h')) return false;
      gIndices[i-2] = _HorP == 'p' ? k - No - 1 : k - 1;
    }
    //LOG(1, "FcidumpReader") << name << ":(" << chemistName << "):"
                            //<< line << std::endl;

    // this will be in physics notation
    {// build up gIndexLens (1, N_1, N_1 * N_2, ..., N_1 *...* N_n-1)
      // copy lens into tempLens
      std::vector<int> tempLens(lens);
      // remove the last lens, the last element
      tempLens.pop_back();
      // insert at the beginning a 1
      tempLens.insert(tempLens.begin(), 1);
      // create the list with the partial products of the elements
      // in tempLens, and store them in gIndexLens.begin()
      std::partial_sum(tempLens.begin(), tempLens.end(),
                       gIndexLens.begin(), std::multiplies<int>());
    }
    // gIndices was read in chemist notation since it comes from a chemistName
    // so we have to change it back to physics notation to store the index
    // correctly
    if (gIndices.size() == 4) gIndices = permuteIndices(gIndices, 1, 2);
    values.push_back(std::atof(std::string{matches[1]}.c_str()));
    indices.push_back(
      std::inner_product(
        gIndexLens.begin(), gIndexLens.end(), gIndices.begin(), 0));
    return true;
  }

  CTF::Tensor<double>* allocateTensor() {
    auto t(new CTF::Tensor<double>(lens.size(),
                                   lens.data(),
                                   syms.data(),
                                   *Cc4s::world));
    t->write(indices.size(), indices.data(), values.data());
    return t;
  }

};

void FcidumpReader::run() {
  const auto filePath(getTextArgument("file", "FCIDUMP"));
  // override the header of the fcidump
  const int nelec(getIntegerArgument("nelec", -1));
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
    {"hhhh", {NO,NO,NO,NO}, "ijkl"},
    {"hhhp", {NO,NO,NO,NV}, "ijka"},
    {"hhph", {NO,NO,NV,NO}, "ijak"},
    {"hhpp", {NO,NO,NV,NV}, "ijab"},
    {"hphh", {NO,NV,NO,NO}, "iajk"},
    {"hphp", {NO,NV,NO,NV}, "iajb"},
    {"hpph", {NO,NV,NV,NO}, "iabj"},
    {"hppp", {NO,NV,NV,NV}, "iabc"},
    {"phhh", {NV,NO,NO,NO}, "aijk"},
    {"phhp", {NV,NO,NO,NV}, "aijb"},
    {"phph", {NV,NO,NV,NO}, "aibj"},
    {"phpp", {NV,NO,NV,NV}, "aibc"},
    {"pphh", {NV,NV,NO,NO}, "abij"},
    {"pphp", {NV,NV,NO,NV}, "abic"},
    {"ppph", {NV,NV,NV,NO}, "abci"},
    {"pppp", {NV,NV,NV,NV}, "abcd"},
  });

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

  for (auto& parser: integralParsers) {
    const auto& name(parser.name);
    if (isArgumentGiven(name)) {
      LOG(0, "FcidumpReader") << "Exporting: " << parser.name << std::endl;
      allocatedTensorArgument<double>(name, parser.allocateTensor());
    }
  }

  std::vector<int> vv({Nv, Nv}), oo({No, No}), ov({No, Nv}), vo({Nv, No});
  std::vector<int> syms({NS, NS});

  // ij: HH HHHH
  auto fij(new CTF::Tensor<double>(2, oo.data(), syms.data(), *Cc4s::world));
  auto hh(getTensorArgument<double>("hh"));
  auto hhhh(getTensorArgument<double>("hhhh"));
  (*fij)["ij"]  = (*hh)["ij"];
  (*fij)["ij"] += (+2.0)*(*hhhh)["ikjk"];
  (*fij)["ij"] += (-1.0)*(*hhhh)["ikkj"];
  allocatedTensorArgument<double>("HHFockMatrix", fij);

  // ab: PP PHPH PHHP
  auto fab(new CTF::Tensor<double>(2, vv.data(), syms.data(), *Cc4s::world));
  auto pp(getTensorArgument<double>("pp"));
  auto phph(getTensorArgument<double>("phph"));
  auto hphp(getTensorArgument<double>("hphp"));
  auto phhp(getTensorArgument<double>("phhp"));
  auto hhpp(getTensorArgument<double>("hhpp"));
  (*fab)["ab"] = (*pp)["ab"];
  (*fab)["ab"] += (+2.0)*(*phph)["akbk"];
  //(*fab)["ab"] += (+2.0)*(*hphp)["kakb"];
  //(*fab)["ab"] += (+2.0)*(*hhpp)["kkab"];
  (*fab)["ab"] += (-1.0)*(*phhp)["akkb"];
  allocatedTensorArgument<double>("PPFockMatrix", fab);

  // ai: PH PHHH
  auto fai(new CTF::Tensor<double>(2, vo.data(), syms.data(), *Cc4s::world));
  auto ph(getTensorArgument<double>("ph"));
  auto phhh(getTensorArgument<double>("phhh"));
  (*fai)["ai"] = (*ph)["ai"];
  (*fai)["ai"] += (+2.0)*(*phhh)["akik"];
  (*fai)["ai"] += (-1.0)*(*phhh)["akki"];
  allocatedTensorArgument<double>("PHFockMatrix", fai);

  // ia: HP HHPH HHHP
  auto fia(new CTF::Tensor<double>(2, ov.data(), syms.data(), *Cc4s::world));
  (*fai)["ai"] = (*fia)["ia"];
  allocatedTensorArgument<double>("HPFockMatrix", fia);

  auto epsi(new CTF::Tensor<double>(1, oo.data(), syms.data(), *Cc4s::world));
  (*epsi)["i"] = (*fij)["ii"];
  allocatedTensorArgument<double>("HoleEigenEnergies", epsi);

  auto epsa(new CTF::Tensor<double>(1, vv.data(), syms.data(), *Cc4s::world));
  (*epsa)["a"] = (*fab)["aa"];
  allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);

  //phph->print();
  //phhp->print();

  CTF::Scalar<double> energy;
  energy[""]  = (+2.0) * (*hh)["ii"];
  energy[""] += (+2.0) * (*hhhh)["ikik"];
  energy[""] += (-1.0) * (*hhhh)["ikki"];
  const double dEnergy(energy.get_val());

  LOG(0, "FcidumpReader") << "hf energy: " << dEnergy << std::endl;

}
