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

}

#endif
