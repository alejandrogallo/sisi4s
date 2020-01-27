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
    {"HHHHCoulombIntegrals", {NO,NO,NO,NO}, "ijkl"},
    {"HHHPCoulombIntegrals", {NO,NO,NO,NV}, "ijka"},
    {"HHPHCoulombIntegrals", {NO,NO,NV,NO}, "ijak"},
    {"HHPPCoulombIntegrals", {NO,NO,NV,NV}, "ijab"},
    {"HPHHCoulombIntegrals", {NO,NV,NO,NO}, "iajk"},
    {"HPHPCoulombIntegrals", {NO,NV,NO,NV}, "iajb"},
    {"HPPHCoulombIntegrals", {NO,NV,NV,NO}, "iabj"},
    {"HPPPCoulombIntegrals", {NO,NV,NV,NV}, "iabc"},
    {"PHHHCoulombIntegrals", {NV,NO,NO,NO}, "aijk"},
    {"PHHPCoulombIntegrals", {NV,NO,NO,NV}, "aijb"},
    {"PHPHCoulombIntegrals", {NV,NO,NV,NO}, "aibj"},
    {"PHPPCoulombIntegrals", {NV,NO,NV,NV}, "aibc"},
    {"PPHHCoulombIntegrals", {NV,NV,NO,NO}, "abij"},
    {"PPHPCoulombIntegrals", {NV,NV,NO,NV}, "abic"},
    {"PPPHCoulombIntegrals", {NV,NV,NV,NO}, "abci"},
    {"PPPPCoulombIntegrals", {NV,NV,NV,NV}, "abcd"},
  });
  std::set<std::string> antiSet;
  for (const auto &integral : infos) {
    if ( ! isArgumentGiven(integral.name) ) continue;
    for (const auto& antigral: integral.getAntisymmetrizers()) {
      if (isArgumentGiven(antigral.name) && !antiSet.count(antigral.name)) {
        // we can use antigral since it has been computed and it has not
        // been initialized
        LOG(1, "CoulombIntegralsFromGaussian")
          << "Anti symmetrizing " << integral.name
          << " with " << antigral.name << std::endl;
        LOG(1, "CoulombIntegralsFromGaussian")
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
      LOG(1, "Integrals")
        << "error: " 
        << integral.name << " could not be antisymmetrized"
        << std::endl;
      throw new EXCEPTION("Antisymmetrization error");
    }
  }


}
