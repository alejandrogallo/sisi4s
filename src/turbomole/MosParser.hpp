#ifndef _TURBOMOLE_MOVECESP_PARSER
#define _TURBOMOLE_MOVECESP_PARSER

#include <regex>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <assert.h>
#include <util/Log.hpp>
#include <util/Parsing.hpp>

using namespace pars;
namespace tmole {

struct MosParser {
  size_t Np;
  std::vector<size_t> Nsaos;
  std::vector<double> occupations, eigenvalues, mos;
  std::smatch match;
  const pars::Regex sep       = blank + oneOrMore
                  , number    = digit + oneOrMore
                  , twodigits = digit + digit
                  , expnumber = anyChar + "." + digit + oneOrMore +
                                "D" + anyChar + twodigits.s
                  , symm      = alnum + oneOrMore
                  , comment   = bof + "#" + print + anyOf
                  , unper     = bof + anyChar + "scfmo" + print + anyOf
                  , eigenval  = bof + sep.s + capture(number.s) + sep.s +
                                capture(symm.s) + sep.s + "eigenvalue=" +
                                capture(expnumber.s) + sep.s + "nsaos=" +
                                capture(number.s) + print + anyOf
                  , mosline1   = capture(expnumber.s)
                  , mosline2   = capture(expnumber.s) + capture(expnumber.s)
                  , mosline3   = capture(expnumber.s) + capture(expnumber.s)
                               + capture(expnumber.s)
                  , mosline4   = capture(expnumber.s) + capture(expnumber.s)
                               + capture(expnumber.s) + capture(expnumber.s)
                  ;

  MosParser(const std::string &fileName) {
    std::fstream file(fileName.c_str(), std::ios::binary | std::ios::in);
    std::vector<double> mos_buf;
    std::string s, e;
    size_t buf_size;
    if (!file.good())    throw "File " + fileName + " not found";
    if (!file.is_open()) throw "File IO error: " + fileName;
    std::string line;
    while (std::getline(file, line)) {

      // empty lines ro comments
      if ( std::regex_match(line, match, comment.r)
        || std::regex_match(line, match, unper.r)
         ) continue;

      // eigenvalues
      if (std::regex_match(line, match, eigenval.r)) {
        e = std::string(match[3]);
        mos_buf.resize(0);
        buf_size = atoi(std::string(match[4]).c_str());
        Nsaos.push_back(buf_size);
      }

      if (std::regex_match(line, match, mosline1.r))
      for (int i(1); i < match.size(); i++) {
        s = std::string(match[i]);
        std::replace( s.begin(), s.end(), 'D', 'e');
        mos_buf.push_back(std::atof(s.c_str()));
      }

      if ( std::regex_match(line, match, mosline2.r))
      for (int i(1); i < match.size(); i++) {
        s = std::string(match[i]);
        std::replace( s.begin(), s.end(), 'D', 'e');
        mos_buf.push_back(std::atof(s.c_str()));
      }

      if ( std::regex_match(line, match, mosline3.r))
      for (int i(1); i < match.size(); i++) {
        s = std::string(match[i]);
        std::replace( s.begin(), s.end(), 'D', 'e');
        mos_buf.push_back(std::atof(s.c_str()));
      }

      if ( std::regex_match(line, match, mosline4.r))
      for (int i(1); i < match.size(); i++) {
        s = std::string(match[i]);
        std::replace( s.begin(), s.end(), 'D', 'e');
        mos_buf.push_back(std::atof(s.c_str()));
      }

      if (mos_buf.size() > 0 && mos_buf.size() == buf_size) {
        std::replace(e.begin(), e.end(), 'D', 'e');
        eigenvalues.push_back(std::atof(e.c_str()));
        for (int64_t i(0); i < mos_buf.size(); i++){
           mos.push_back(mos_buf[i]);
        }
        mos_buf.resize(0);
      }
    }

    Np = eigenvalues.size();

    file.close();

  }
};


}

#endif
