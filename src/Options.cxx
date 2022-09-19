/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <Options.hpp>
#include <string>
#include <sstream>

using namespace sisi4s;

Options::Options(int argc_, char **argv_)
  : app{"SiSi4S: Coupled Cluster For Solids"}
  , logLevel(Options::DEFAULT_LOG_LEVEL)
  , inFile("")
  , logFile("sisi4s.log")
  , yamlOutFile("sisi4s.out.yaml")
  , argc(argc_)
  , argv(argv_)
  , cc4s(false)
  , listAlgorithms(false)
  , dryRun(false)
  {

    app.add_option("-i,--in", inFile, "Input file path")
      ->check(CLI::ExistingFile)
      ->required(true);

    app.add_option("-o,--out", yamlOutFile, "Output yaml file")
      ->default_val(yamlOutFile);

    app.add_option("-l,--log", logFile, "Output log file")
      ->default_val(logFile);

    app.add_flag("--dry", dryRun, "Do a dry run pass")
      ->default_val(dryRun);

    app.add_option("--log-level", logLevel, "Log level")
      ->default_val(logLevel);

    app.add_flag("--cc4s", cc4s,
                 "Interpret the input file in the old cc4s DSL")
      ->default_val(cc4s);

    app.add_flag("--list-algorithms,--list", listAlgorithms,
                 "List registered algorithms")
      ->default_val(listAlgorithms);

}

int Options::parse() {
  try {
    (app).parse((argc), (argv));
  } catch(const CLI::ParseError &e) {
    (app).exit(e);
    std::exit(1);
  }
  if (app.get_help_ptr()->count()) std::exit(0);
  return 0;
}
