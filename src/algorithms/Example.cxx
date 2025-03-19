#include <Step.hpp>
using namespace sisi4s;

#include <SpecChain.hpp>

DEFSPEC(
    Example,
    SPEC_IN(
        {"name",
         SPEC_CHAIN("Example of a input that can be a variable name or a value",
                    SPEC_VARIN("doc", std::string),
                    SPEC_VALUE("doc", std::string))
             ->require()},
        {"suffix", SPEC_VALUE_DEF("Suffix to make", std::string, "etc")},
        {"range", SPEC_RANGE("A range example", int, 1, 13)},
        {"vector",
         SPEC_VALUE_DEF("An example for a vector of numbers",
                        std::vector<int>,
                        {1, 2, 3})},
        {"age", SPEC_POSITIVE("Positive number", int)},
        {"chain",
         SPEC_CHAIN(
             "An example of a chained spec, it can be double, int or string",
             SPEC_NEGATIVE("Negative double number", double),
             SPEC_VALUE_DEF("", int, 42),
             SPEC_VALUE_DEF("", std::string, "Whatever"))}),
    SPEC_OUT({"name", SPEC_VAROUT("doc", std::string)}));

DEFSTEP(Example) {
  const auto name = in.get<std::string>("name"),
             suffix = in.get<std::string>("suffix"), out_name = name + suffix;

  if (in.present("chain")) {
    if (in.is_of_type<double>("chain")) {
      std::cout << "chain <double>: " << in.get<double>("chain") << std::endl;
    } else if (in.is_of_type<int>("chain")) {
      std::cout << "chain <int>: " << in.get<int>("chain") << std::endl;
    } else if (in.is_of_type<std::string>("chain")) {
      std::cout << "chain <str>: " << in.get<std::string>("chain") << std::endl;
    }
  }
  std::cout << "Name : " << name << std::endl;
  std::cout << "OutName : " << out_name << std::endl;

  if (out.present("name")) out.set<std::string>("name", out_name);
}
