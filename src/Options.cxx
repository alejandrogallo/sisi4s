#include <Options.hpp>
#include <string>
#include <sstream>

using namespace sisi4s;

Options::Options(int argc_, char **argv_)
    : app{"SiSi4S: Coupled Cluster For Solids"}
    , log_level(Options::DEFAULT_LOG_LEVEL)
    , in_file("")
    , log_file("sisi4s.log")
    , yaml_out_file("sisi4s.out.yaml")
    , lisp_file("sisi.yaml")
    , name("")
    , argc(argc_)
    , argv(argv_)
    , list_algorithms_p(false)
    , dryRun(false)
    , algo_specs({}) {

  app.add_option("-i,--in", in_file, "Input file path")
      ->check(CLI::ExistingFile);

  app.add_option("-n,--name", name, "Name of the calculation");

  app.add_option("-o,--out", yaml_out_file, "Output yaml file")
      ->default_val(yaml_out_file);

  app.add_option("-l,--log", log_file, "Output log file")
      ->default_val(log_file);

  app.add_flag("--dry", dryRun, "Do a dry run pass")->default_val(dryRun);
  app.add_flag("-f,--force",
               force,
               "Run calculation even if there are warnings in the input");
  app.add_flag("-c,--only-check-input",
               only_check_input,
               "Check only the validity of the input file");

  app.add_option("--log-level", log_level, "Log level")->default_val(log_level);

#if defined(HAVE_LISP)
  app.add_option("--load", lisp_file, "Lisp Input file path")
      ->check(CLI::ExistingFile);
#endif /* defined(HAVE_LISP) */

  app.add_flag("--list-algorithms,--list",
               list_algorithms_p,
               "List registered algorithms")
      ->default_val(list_algorithms_p);

  app.add_option("--steps",
                 algo_specs,
                 "The name(s) of the step you want to consult");
}

int Options::parse() {
  try {
    (app).parse((argc), (argv));
  } catch (const CLI::ParseError &e) {
    (app).exit(e);
    std::exit(1);
  }
  if (app.get_help_ptr()->count()) std::exit(0);
  if (name.size()) {
    log_file = name + ".log";
    yaml_out_file = name + ".out.yaml";
  }
  return 0;
}
