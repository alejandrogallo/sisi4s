#ifndef _NUCLEAR_REPULSION_ENERGY_DEFINED
#define _NUCLEAR_REPULSION_ENERGY_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

class NuclearRepulsionEnergy : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(NuclearRepulsionEnergy);
  NuclearRepulsionEnergy(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~NuclearRepulsionEnergy() {}
  virtual void run();
};
} // namespace sisi4s

#endif
