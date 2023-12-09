#include <fstream>
#include <locale>

#include <AlgorithmInputSpec.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Exception.hpp>
#include <yaml-cpp/yaml.h>
#include <NewData.hpp>

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

    if (node["disable"] && node["disable"].as<bool>()) continue;
    if (node["enable"] && !node["enable"].as<bool>()) continue;

    std::string name = node["name"].as<std::string>();
    Arguments in, out;
    const auto spec = AlgorithmInputSpec::map[name];

    struct Phase {
      char const *section;
      AlgorithmInputSpec::Spec const &spec;
      Arguments &args;
    };
    std::vector<Phase> phases = {{"in", spec.in, in}, {"out", spec.out, out}};

    std::cout << "\n\nCalling spec for " << name << ">> " << std::endl;

    for (auto &phase : phases) {
      for (auto const &pair : phase.spec) {
        auto spec = pair.second;
        const std::string key_name = Arguments::normalize_name(pair.first);
        YAML::Node input_settings = node[phase.section], yaml_node;
        bool found_node_p = false;
        std::string provided_key;
        for (YAML::Node::const_iterator it = input_settings.begin();
             it != input_settings.end();
             it++) {
          provided_key = it->first.as<std::string>();
          if (Arguments::normalize_name(provided_key) == key_name) {
            yaml_node = it->second;
            found_node_p = true;
          }
        }
        /* const std::string provided_key = found_pair->second; */
        if (found_node_p) {
          const std::string yaml_string = yaml_node.Scalar();
          std::cout << "\t> " << name << "." << phase.section << "."
                    << pair.first << std::endl;
          spec->parse(yaml_string);
          if (!spec->validate()) {
            std::cout << "\t\tnot valid !  " << spec->validate() << std::endl;
            std::cout << "\t\t: Doc: " << spec->doc << std::endl;
            const auto warnings = spec->warnings(yaml_string);
            std::cout << "\t\t✗ ERROR checking spec (" << warnings.size()
                      << " warnings)" << std::endl;
            for (auto const &warning : warnings) {
              std::cout << "\t\t- " << warning << std::endl;
              std::exit(1);
            }
          }
          const std::string db_index = spec->commit();
          phase.args.push(key_name, db_index);
        } else {
          if (spec->has_default) {
            std::cout << "\t> (def) " << name << "." << phase.section << "."
                      << pair.first << std::endl;
            phase.args.push(key_name, spec->commit());
          }
        }
        // check foreign keys
      }
    }

    Algorithm *algorithm(AlgorithmFactory::create(name, in, out));
    if (node["note"]) { algorithm->note = node["note"].as<std::string>(); }
    if (node["fallible"]) { algorithm->fallible = node["fallible"].as<bool>(); }
    algorithms.push_back(algorithm);
  }
  return algorithms;
}
