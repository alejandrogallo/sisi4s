#include <algorithms/TensorAntisymmetrizer.hpp>
#include <string>
#include <vector>
#include <math/MathFunctions.hpp>
#include <algorithm>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
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
  std::map<std::string, PTR(CTF::Tensor<double>)> integralCopies;
  for (const auto &integral : infos) {
  if ( isArgumentGiven(integral.name) ) {
    LOG(1, "TensorAntisymmetrizer") << "Copying " << integral.name
      << std::endl;
    integralCopies[integral.name] = NEW(
      CTF::Tensor<double>,
      *getTensorArgument<double>(integral.name));
  }
  }

#ifdef DEBUG
  for (auto& copy: integralCopies) {
    double norm(frobeniusNorm(*copy.second));
    LOG(1, "TensorAntisymmetrizerNorm") << "|" << copy.first << "| = "
      << norm << std::endl;
  }
#endif

  for (const auto &integral : infos) {
    bool antisymmetrized(false);
    if ( ! isArgumentGiven(integral.name) ) continue;

    for (const auto& antigral: integral.getAntisymmetrizers()) {
      if (isArgumentGiven(antigral.name)) {
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
        //auto aantiCtf(getTensorArgument<double>(antigral.name));
        auto antiCtf(integralCopies[antigral.name]);
        (*inteCtf)[integral.ids.data()] -= (*antiCtf)[antigral.ids.data()];
        antisymmetrized = true;
        // go to the next integral
        break;
      }
    }

    // if integral.name is not antisymmetrized, then something went
    // wrong, we don't have enough integrals to do this, panic!
    if (!antisymmetrized) {
      LOG(1, "TensorAntisymmetrizer")
        << "error: " 
        << integral.name << " could not be antisymmetrized"
        << std::endl;
      LOG(1, "TensorAntisymmetrizer")
        << "error: possibilities: " << std::endl;
      for (const auto& antigral: integral.getAntisymmetrizers()) {
        LOG(1, "TensorAntisymmetrizer") << "> " << antigral.name << std::endl;
      }
      throw new EXCEPTION("Antisymmetrization error");
    }
  }

#ifdef DEBUG
  for (auto& copy: integralCopies) {
    double norm(frobeniusNorm(*copy.second));
    LOG(1, "TensorAntisymmetrizerNorm") << "|" << copy.first << "| = "
      << norm << std::endl;
  }
#endif

}
