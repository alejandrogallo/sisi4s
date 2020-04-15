#include <vector>
#include <algorithms/OneBodyRotation.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <ctf.hpp>
#include <util/Emitter.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(OneBodyRotation);
#define LOGGER(_l) LOG(_l, "OneBodyRotation")

void OneBodyRotation::run() {

  checkArgumentsOrDie( { "OrbitalCoefficients"
                       , "Data", "No"
                       , "Out", "hh", "pp", "hp", "ph"
                       } );

  auto C(getTensorArgument("OrbitalCoefficients"));
  auto I(getTensorArgument("Data"));
  auto O(new CTF::Tensor<double>(*I));
  const int No(getIntegerArgument("No", 0)), Np(C->lens[0]);

  LOGGER(0) << "No: " << No << std::endl;
  LOGGER(0) << "Np: " << Np << std::endl;

  LOGGER(0) << "Rotating" << std::endl;
  (*O)["pq"] = (*C)["Pp"]
             * (*C)["Qq"]
             * (*I)["PQ"]
             ;

  struct _info { const char *name; const int begin[2], end[2]; };
  const std::vector<_info> infos = { {  "hh", {0,   0}, {No, No} }
                                   , {  "pp", {No, No}, {Np, Np} }
                                   , {  "hp", {0,  No}, {No, Np} }
                                   , {  "ph", {No,  0}, {Np, No} }
                                   , { "Out", {0,   0}, {Np, Np} }
                                   };

  for (const auto& i: infos) {
    if ( ! isArgumentGiven(i.name) ) continue;
    auto tensor(new CTF::Tensor<double>(O->slice(i.begin, i.end)));
    LOGGER(0) << i.name << ": "
              << "{" << i.begin[0] << "," << i.begin[1] << "}"
              << " --> "
              << "{" << i.end[0]   << "," << i.end[1]   << "}"
              << std::endl;
    allocatedTensorArgument<double>(i.name, tensor);
  }

}
