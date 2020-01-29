#include <algorithms/TensorAntisymmetrizer.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <map>
#include <set>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorAntisymmetrizer);

void TensorAntisymmetrizer::run() {
  const std::vector<IntegralInfo> infos({
    {"HPHPCoulombIntegrals", {NO,NV,NO,NV}, "iajb"},
    {"PHHPCoulombIntegrals", {NV,NO,NO,NV}, "aijb"},
    {"PHPHCoulombIntegrals", {NV,NO,NV,NO}, "aibj"},
    {"HPPHCoulombIntegrals", {NO,NV,NV,NO}, "iabj"},
    // trivially self antisymmetrizable
    {"HHPHCoulombIntegrals", {NO,NO,NV,NO}, "ijak"},
    {"HHPPCoulombIntegrals", {NO,NO,NV,NV}, "ijab"},
    {"PPHPCoulombIntegrals", {NV,NV,NO,NV}, "abic"},
    {"PPHHCoulombIntegrals", {NV,NV,NO,NO}, "abij"},
    {"HHHHCoulombIntegrals", {NO,NO,NO,NO}, "ijkl"},
    {"HHHPCoulombIntegrals", {NO,NO,NO,NV}, "ijka"},
    {"PPPPCoulombIntegrals", {NV,NV,NV,NV}, "abcd"},
    {"PPPHCoulombIntegrals", {NV,NV,NV,NO}, "abci"},
    {"HPHHCoulombIntegrals", {NO,NV,NO,NO}, "iajk"},
    {"PHHHCoulombIntegrals", {NV,NO,NO,NO}, "aijk"},
    {"HPPPCoulombIntegrals", {NO,NV,NV,NV}, "iabc"},
    {"PHPPCoulombIntegrals", {NV,NO,NV,NV}, "aibc"},
  });
  std::set<std::string> antiSet;
  for (const auto &integral : infos) {
    if ( ! isArgumentGiven(integral.name) ) continue;
    for (const auto& antigral: integral.getAntisymmetrizers()) {
      if (isArgumentGiven(antigral.name) && !antiSet.count(antigral.name)) {
        // we can use antigral since it has been computed and it has not
        // been initialized
        LOG(1, "TensorAntisymmetrizer")
          << integral.name
          << " from " << antigral.name << std::endl;
        LOG(1, "TensorAntisymmetrizer") << "  "
          << integral.name << "[" << integral.ids << "] -= "
          << antigral.name << "[" << antigral.ids << "]"
          << std::endl;
        auto inteCtf(getTensorArgument<double>(integral.name));
        auto antiCtf(getTensorArgument<double>(antigral.name));
        (*inteCtf)[integral.ids.data()] -= (*antiCtf)[antigral.ids.data()];
        antiSet.insert(integral.name);
        break;
      }

    }
    // if integral.name is not antisymmetrized, then something went
    // wrong, we don't have enough integrals to do this, panic!
    if (!antiSet.count(integral.name)) {
      LOG(1, "TensorAntisymmetrizer")
        << "error: " 
        << integral.name << " could not be antisymmetrized"
        << std::endl;
      LOG(1, "TensorAntisymmetrizer")
        << "error: I tried the following: " << std::endl;
      for (const auto& antigral: integral.getAntisymmetrizers()) {
        LOG(1, "TensorAntisymmetrizer") << "> " << antigral.name << std::endl;
      }
      throw new EXCEPTION("Antisymmetrization error");
    }
  }


}
