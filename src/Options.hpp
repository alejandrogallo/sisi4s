#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

#include <string>
#include <vector>
#include <vendor/CLI11.hpp>
#include "config.h"

namespace sisi4s {
class Options {
public:
  CLI::App app;
  int log_level;
  std::string in_file, log_file, yaml_out_file, lisp_file, name;
  int argc;
  char **argv;
  bool list_algorithms_p, dryRun;
  std::vector<std::string> algo_specs;

  static const int DEFAULT_LOG_LEVEL = 1;

  Options(int argc, char **argv);
  int parse();
};
} // namespace sisi4s

#endif
