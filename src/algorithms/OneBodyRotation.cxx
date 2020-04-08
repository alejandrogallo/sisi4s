#include <string>
#include <vector>
#include <algorithm>
#include <algorithms/OneBodyRotation.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>
#include <util/Emitter.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(OneBodyRotation);
#define LOGGER(_l) LOG(_l, "OneBodyRotation")

void OneBodyRotation::run() {

  checkArgumentsOrDie( { "OrbitalCoefficients"
                       , "Data"
                       , "Out"
                       } );

  auto C(getTensorArgument("OrbitalCoefficients"));
  auto I(getTensorArgument("Data"));
  auto O(new CTF::Tensor<double>(*I));

  LOGGER(0) << "Computing... " << std::endl;
  (*O)["pq"] = (*C)["Pp"]
             * (*C)["Qq"]
             * (*I)["PQ"]
             ;

  allocatedTensorArgument<double>("Out", O);

}
