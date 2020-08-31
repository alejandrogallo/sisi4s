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
// p|q
// r s
template <typename A>
A upTr(const A& a) { return permuteIndices(a, 0, 1); }
// p q
// r|s
template <typename A>
A downTr(const A& a) { return permuteIndices(a, 2, 3); }

// struct used to store information about the names, indices ranges
// and string indices for ctf. It computes also IntegralInfos for
// antisymmetry options
struct IntegralInfo {
  std::string name;
  std::array<Index, 4> indices;
  std::string ids;
  IntegralInfo(std::string n, std::array<Index, 4> i, std::string is):
    name(n), indices(i), ids(is){}

  std::vector<IntegralInfo> getAntisymmetrizers() const {
    std::vector<IntegralInfo> result;
    std::vector<IntegralInfo> targets;
    std::set<std::string> names;
    // build up targets
    targets.push_back({upTr(name), upTr(indices), upTr(ids)});
    targets.push_back({downTr(name), downTr(indices), downTr(ids)});
    for (auto& t: targets) {
      std::vector<IntegralInfo> equivalents;
      // set up the equivalent tensors for the current target
      equivalents.push_back({vSym(t.name),  vSym(t.indices),  vSym(t.ids)});
      equivalents.push_back({hSym(t.name),  hSym(t.indices),  hSym(t.ids)});
      equivalents.push_back({vlSym(t.name), vlSym(t.indices), vlSym(t.ids)});
      for (auto& eq: equivalents) {
        if (names.count(eq.name) >= 1) continue;
        names.insert(eq.name);
        result.push_back({eq.name, eq.indices, eq.ids});
      }
    }
    return result;
  }

};

};

#endif
