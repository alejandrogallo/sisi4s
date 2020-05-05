#include <util/Parsing.hpp>
#include <vector>
#include <string>

namespace pars {

std::string oneOf(const std::vector<std::string> &v) {
  return std::accumulate
    ( v.begin()
    , v.end()
    , std::string("")
    , std::function<Str(Str, Str)>(
        [](const Str &a, const Str &b){ return a.size() ? a + orOf + b : b; }
        )
    )
    ;
}

template <>
std::vector<std::string> parseVector(const std::string &l) {
  const Regex base = capture("[^,]" + oneOrMore);
  std::smatch match;
  std::string lcopy(l);
  std::vector<std::string> result;
  while(std::regex_search(lcopy, match, base.r)) {
    result.push_back(match[1].str());
    lcopy = match.suffix();
  }
  return result;
}

template<>
std::vector<double> parseVector<double>(const std::string &l) {
  std::vector<double> res;
  for (const auto &s: parseVector<std::string>(l))
    res.push_back(std::atof(s.c_str()));
  return res;
}

template<>
std::vector<int> parseVector<int>(const std::string &l) {
  std::vector<int> res;
  for (const auto &s: parseVector<std::string>(l))
    res.push_back(std::atoi(s.c_str()));
  return res;
}

}
