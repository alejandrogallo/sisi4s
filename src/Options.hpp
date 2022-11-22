/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

#include <string>
#include <vendor/CLI11.hpp>

namespace sisi4s {
struct Options {

  CLI::App app;
  int logLevel;
  std::string inFile, logFile, yamlOutFile;
  int argc;
  char **argv;
  bool cc4s, listAlgorithms, dryRun;

  static const int DEFAULT_LOG_LEVEL = 1;

  Options(int argc, char **argv);
  int parse();
};
} // namespace sisi4s

#endif
