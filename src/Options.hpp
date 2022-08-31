/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

#include <string>
#include <extern/CLI11.hpp>

namespace cc4s {
  struct Options {

    int logLevel;
    std::string inFile, logFile, yamlOutFile;
    bool dryRun;
    int argc;
    char** argv;

    const int DEFAULT_LOG_LEVEL = 1;
    bool hummel;
    CLI::App app;

    Options(int argc, char **argv);
    int parse();
  };
}

#endif

