#include "Emitter.hpp"

#include <fstream>
#include <string>

using namespace sisi4s;

int Emitter::rank(-1);
std::string Emitter::fileName("sisi4s.yaml");
PTR(std::ofstream) Emitter::yamlFile;
PTR(YAML::Emitter) Emitter::yamlEmitter;

void Emitter::setRank(int const rank_) { rank = rank_; }

int Emitter::getRank() { return rank; }

void Emitter::setFileName(const std::string &name) { fileName = name; }

YAML::Emitter &Emitter::getEmitter() {
  if (!yamlFile) {
    yamlFile =
        NEW(std::ofstream, fileName, std::ofstream::out | std::ofstream::trunc);
  }
  if (!yamlEmitter) { yamlEmitter = NEW(YAML::Emitter, *yamlFile); }
  return *yamlEmitter;
}

void Emitter::flush() { std::flush(*yamlFile); }
