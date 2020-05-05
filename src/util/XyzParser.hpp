#ifndef _XYZ_PARSER_
#define _XYZ_PARSER_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <util/Parsing.hpp>
#include <regex>


struct Atom {
  const std::string symbol;
  const struct { double x, y, z; } position;
};

namespace pars {

struct XyzParser {

  const Regex sep = blank + oneOrMore // any number > 1 of spaces or tabs
            , atom = upper + lower + optional // atom is Upper + lower?
            , xyz_line = bof + blank + anyOf // spaces at the start allowed
                       + capture(atom.s) + sep.s // capture atom symbol
                       + capture(realNumber) + sep.s
                       + capture(realNumber) + sep.s
                       + capture(realNumber)
                       + blank + anyOf + eof
            ;

  Atom parseLine(const std::string &line) {
    std::smatch m;
    std::regex_match(line, m, xyz_line.r);
    if (m.size() == 0) { throw "Line " + line + " is malformed"; }
    return { std::string(m[1])
           , { std::atof(std::string(m[2]).c_str())
             , std::atof(std::string(m[3]).c_str())
             , std::atof(std::string(m[4]).c_str())
             }
           };
  }

  std::vector<Atom> parseFile(const std::string &fileName) {
    std::fstream f(fileName);
    std::string line;
    std::vector<Atom> atoms;
    if (!f.is_open()) throw "File IO error: " + fileName;

    std::getline(f, line); // Number of atoms
    const int natoms = std::atoi(line.c_str());

    std::getline(f, line); // Title of xyz file

    // parse lines
    for (int i(0); i<natoms; i++){
      std::getline(f, line);
      atoms.push_back(parseLine(line));
    }

    return atoms;
  }

};

}

#endif
