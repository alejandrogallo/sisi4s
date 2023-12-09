#include <Parser.hpp>

#include <algorithms/Algorithm.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <locale>
#include <yaml-cpp/yaml.h>

using namespace sisi4s;

///////////////////////////////////////////////////////////////////////////////
// YAML Parser ////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

InputFileParser<InputFileFormat::YAML>::InputFileParser(
    std::string const &fileName_)
    : fileName(fileName_) {}

InputFileParser<InputFileFormat::YAML>::~InputFileParser() {}

std::vector<Algorithm *> InputFileParser<InputFileFormat::YAML>::parse() {
  const YAML::Node nodes = YAML::LoadFile(fileName);
  std::vector<Algorithm *> algorithms;

  for (const YAML::Node &node : nodes) {
    std::string name = node["name"].as<std::string>();
    std::vector<Argument> arguments;

    if (node["disable"] && node["disable"].as<bool>()) continue;
    if (node["enable"] && !node["enable"].as<bool>()) continue;

    for (auto const &_inout : {"in", "out"}) {
      YAML::Node const inout = node[_inout];
      for (YAML::const_iterator it = inout.begin(); it != inout.end(); ++it) {
        std::string key = it->first.as<std::string>(), valueName;
        try {
          int value = it->second.as<int>();
          valueName = (new IntegerData(value))->getName();
        } catch (YAML::TypedBadConversion<int> const &c) {

          try {
            double value = it->second.as<double>();
            valueName = (new RealData(value))->getName();
          } catch (YAML::TypedBadConversion<double> const &c) {

            try {
              std::string value = it->second.as<std::string>();
              if (value.substr(0, 1) == "$") {
                const std::string symbolName = value.substr(1);
                Data *data(Data::get(symbolName));
                valueName =
                    data ? data->getName() : (new Data(symbolName))->getName();
              } else {
                valueName = (new TextData(value))->getName();
              }
            } catch (YAML::TypedBadConversion<std::string> const &c) {
              throw c;
            }
          }
        }

        arguments.push_back(Argument(key, valueName));
      }
    }

    Algorithm *algorithm(AlgorithmFactory::create(name, arguments));
    if (node["note"]) { algorithm->note = node["note"].as<std::string>(); }
    if (node["fallible"]) { algorithm->fallible = node["fallible"].as<bool>(); }
    algorithms.push_back(algorithm);
  }
  return algorithms;
}
