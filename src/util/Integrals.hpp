#ifndef __UTIL_INTEGRALS_DEFINED
#define __UTIL_INTEGRALS_DEFINED
#include <vector>
#include <array>
namespace cc4s {

// general enum for calculating different parts of the integrals
enum Index { NO, NV, NP };

template <typename A>
A permuteIndices(const A& a, size_t i, size_t j) {
  A b(a); b[i] = a[j]; b[j] = a[i]; return b;
}

// struct used to store information about the names, indices ranges
// and string indices for ctf. It computes also IntegralInfos for
// antisymmetry options
struct IntegralInfo {
  const std::string name;
  const std::array<Index, 4> indices;
  const std::string ids;
  IntegralInfo(std::string n, std::array<Index, 4> i, std::string is):
    name(n), indices(i), ids(is){}
  std::vector<IntegralInfo> getAntisymmetrizers() const {
    std::vector<IntegralInfo> result;
    for (auto i: std::vector<int>({0,2})) {
      std::string newName(permuteIndices(name, i, i+1));
      std::array<Index, 4> newIndices(permuteIndices(indices, i, i+1));
      std::string newIds(permuteIndices(ids, i, i+1));
      result.push_back({newName, newIndices, newIds});
    }
    return result;
  }
};

};

#endif
