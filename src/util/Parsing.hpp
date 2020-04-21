#ifndef _ALEJANDROS_PARSING_LIBRAY
#define _ALEJANDROS_PARSING_LIBRAY

#include<regex>
#include<string>
#include<vector>

namespace pars {
  using Str = std::string;
  const
  Str oneOrMore("+")
    , newline("\\n")
    , anyChar(".")
    , tab("\\t")
    , anyOf("*")
    , optional("?")
    , eof("$")
    , orOf("|")
    , bof("^") // begin of line
    , alnum("[[:alnum:]]") // lowercase letters, uppercase letters, and digits
    , alpha("[[:alpha:]]") // lowercase letters and uppercase letters
    , blank("[[:blank:]]") // space or tab
    , cntrl("[[:cntrl:]]") // the file format escape characters
    , digit("[[:digit:]]") // digits
   // lowercase letters, uppercase letters, digits, and punctuation
    , graph("[[:graph:]]")
    , lower("[[:lower:]]") // lowercase letters
    , upper("[[:upper:]]")  // uppercase characters
   // lowercase letters, uppercase letters, digits, punctuation, and space
    , print("[[:print:]]")
    , punct("[[:punct:]]") // punctuation
    , space("[[:space:]]") // space
    , xdigit("[[:xdigit:]]")  // hexadecimal digits (with upper and lower)
    , realNumber("[-eE+\\d.]+") // real number TODO: improve
    ;

  const
  std::function<Str(Str)> capture([](const Str& i) { return "(" + i + ")"; })
                        , group([](const Str& i) { return "(?:" + i + ")"; })
                        ;

  Str oneOf(const std::vector<Str> &v) {
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

  struct Regex {
    const std::string s;
    const std::regex r;
    Regex(const char* s_): s(s_), r(s_) {}
    Regex(const std::string &s_): s(s_), r(s_) {}
  };

  template <typename F> std::vector<F> parseVector(const std::string&);
  template<>
  std::vector<std::string> parseVector(const std::string& l) {
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
  std::vector<double> parseVector(const std::string& l) {
    std::vector<double> res;
    for (const auto &s: parseVector<std::string>(l))
      res.push_back(std::atof(s.c_str()));
    return res;
  }
  template<typename F=int>
  std::vector<F> parseVector(const std::string& l) {
    std::vector<F> res;
    for (const auto &s: parseVector<std::string>(l))
      res.push_back(std::atoi(s.c_str()));
    return res;
  }


}

#endif
