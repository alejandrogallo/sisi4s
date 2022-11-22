#ifndef DEFINEHOLESANDPARTICLES_HPP_
#define DEFINEHOLESANDPARTICLES_HPP_

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
namespace cc4s {

// the header that should be present
// in such a tensor
struct HPHeader {
  double fermiEnergy;
  std::vector<double> energies;
};

} // namespace cc4s

DEFINE_ALGORITHM_HEADER(DefineHolesAndParticles, );

} // namespace sisi4s

#endif
