/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argc_, char **argv_)
  : inFile("cc4s.in")
  , logFile("cc4s.log")
  , yamlOutFile("cc4s.out.yaml")
  , app{"S4CC: Coupled Cluster For Solids"}
  , argc(argc_)
  , argv(argv_)
  , dryRun(false)
  , hummel(false)
  {
  logLevel = DEFAULT_LOG_LEVEL;
  app.add_option("-i,--in", inFile, "Input file path")
     ->default_val(inFile);
  app.add_option("-o,--out", yamlOutFile, "Output yaml file")
     ->default_val(yamlOutFile);
  app.add_option("-l,--log", logFile, "Output log file")
     ->default_val(logFile);
  app.add_option("--dry", dryRun, "Do a dry run pass")
     ->default_val(dryRun);
  app.add_option("--log-level", logLevel, "Log level")
     ->default_val(logLevel);
  app.add_option("--hummel", hummel,
                 "Interpret the input file in the old cc4s DSL")
     ->default_val(hummel);
}

int Options::parse() {
  CLI11_PARSE(app, argc, argv);
  return 0;
}

