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

// p | q
// r | s
template <typename A>
A vSym(const A& a) { return permuteIndices(permuteIndices(a, 0, 1), 2, 3); }
// p  q
// -----
// r  s
template <typename A>
A hSym(const A& a) { return permuteIndices(permuteIndices(a, 0, 2), 1, 3); }
// p  q
// -
// r  s
template <typename A>
A vlSym(const A& a) { return permuteIndices(a, 0, 2); }

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
    std::vector<IntegralInfo> equivalents;
    std::set<std::string> names;
    equivalents.push_back({vSym(name), vSym(indices), vSym(ids)});
    equivalents.push_back({hSym(name), hSym(indices), hSym(ids)});
    equivalents.push_back({vlSym(name), vlSym(indices), vlSym(ids)});
    for (auto& eq: equivalents) {
      for (auto i: std::vector<int>({0,2})) {
        std::string newName(permuteIndices(eq.name, i, i+1));
        if (names.count(newName) >= 1) continue;
        names.insert(newName);
        std::array<Index, 4> newIndices(permuteIndices(eq.indices, i, i+1));
        std::string newIds(permuteIndices(eq.ids, i, i+1));
        result.push_back({newName, newIndices, newIds});
      }
    }
    // TODO: sort and unique result
    return result;
  }
  private:

};

};

#endif
