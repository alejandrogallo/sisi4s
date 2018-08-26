#include <algorithms/HartreeFock.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <algorithms/Algorithm.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(HartreeFock);
HartreeFock::HartreeFock(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
HartreeFock::~HartreeFock() {}

void HartreeFock::run() {
  LOG(0, "HartreeFock") << "Calculating hartree fock" << std::endl;
}
