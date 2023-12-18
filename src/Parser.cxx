#include <fstream>
#include <locale>

#include <yaml-cpp/yaml.h>

#include <AlgorithmInputSpec.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <NewData.hpp>

using namespace sisi4s;

///////////////////////////////////////////////////////////////////////////////
// YAML Parser ////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

InputFileParser<InputFileFormat::YAML>::InputFileParser(
    std::string const &fileName_,
    bool e)
    : fileName(fileName_)
    , exit_on_warnings(e) {}

InputFileParser<InputFileFormat::YAML>::~InputFileParser() {}

std::vector<Algorithm *> InputFileParser<InputFileFormat::YAML>::parse() {
  const YAML::Node nodes = YAML::LoadFile(fileName);
  std::vector<Algorithm *> algorithms;
  struct Error {
    std::string step, phase, key;
    int line, column;
    std::vector<std::string> warnings;
  };
  std::vector<Error> errors;

  for (const YAML::Node &node : nodes) {

    if (node["disable"] && node["disable"].as<bool>()) continue;
    if (node["enable"] && !node["enable"].as<bool>()) continue;

    std::string name = node["name"].as<std::string>();
    Arguments in, out;
    const auto spec =
        AlgorithmInputSpec::map[AlgorithmFactory::normalize_name(name)];

    struct Phase {
      char const *section;
      AlgorithmInputSpec::Spec const &spec;
      Arguments &args;
    };
    std::vector<Phase> phases = {{"in", spec.in, in}, {"out", spec.out, out}};

    // std::cout << "\n\nCalling spec for " << name << ">> " << std::endl;

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
            break;
          }
        }
        /* const std::string provided_key = found_pair->second; */
        if (found_node_p) {
          const std::string yaml_string = yaml_node.Scalar();
          // std::cout << "\t> " << name << "." << phase.section << "."
          //           << pair.first << std::endl;
          spec->parse(yaml_string);
          if (!spec->validate()) {
            const auto warnings = spec->warnings(yaml_string);
            // std::cout << "\t\tnot valid !  " << spec->validate() <<
            // std::endl; std::cout << "\t\t: Doc: " << spec->doc << std::endl;
            // std::cout << "\t\t✗ ERROR checking spec (" << warnings.size()
            //           << " warnings)" << std::endl;
            errors.push_back({name,
                              phase.section,
                              provided_key,
                              yaml_node.Mark().line,
                              yaml_node.Mark().column,
                              warnings});
          }
          const std::string db_index = spec->commit();
          phase.args.push(key_name, db_index);
        } else {
          if (spec->has_default) {
            // std::cout << "\t> (def) " << name << "." << phase.section << "."
            //           << pair.first << std::endl;
            phase.args.push(key_name, spec->commit());
          } else if (spec->required) {
            errors.push_back(
                {name,
                 phase.section,
                 pair.first,
                 input_settings.Mark().line,
                 input_settings.Mark().column,
                 {_FORMAT("The key %s is required.", pair.first.c_str())}});
          }
        }
      }

      // Check for foreign keys errors
      YAML::Node input_settings = node[phase.section], yaml_node;
      for (YAML::Node::const_iterator it = input_settings.begin();
           it != input_settings.end();
           it++) {
        const std::string provided_key = it->first.as<std::string>();
        bool found = false;
        for (auto const &sit : phase.spec) {
          const std::string key_name = Arguments::normalize_name(sit.first);
          if (Arguments::normalize_name(provided_key) == key_name) {
            found = true;
            break;
          }
        }
        if (!found) {
          errors.push_back(
              {name,
               phase.section,
               provided_key,
               it->first.Mark().line,
               it->first.Mark().column,
               {_FORMAT("The key %s is not part of the spec and is not "
                        "understood",
                        provided_key.c_str())}});
        }
      }
    }

    Algorithm *algorithm(AlgorithmFactory::create(name, in, out));
    if (node["note"]) { algorithm->note = node["note"].as<std::string>(); }
    if (node["fallible"]) { algorithm->fallible = node["fallible"].as<bool>(); }
    if (algorithm == nullptr) {
      errors.push_back(
          {name,
           "",
           "",
           node.Mark().line,
           node.Mark().column,
           {_FORMAT(
               "The algorithm with name %s is not defined, check the spelling",
               name.c_str())}});
    } else {
      algorithms.push_back(algorithm);
    }
  }

  if (errors.size()) {
    LOG(0, "Parser") << _FORMAT("(%ld) ERRORS Encountered in the input file:\n",
                                errors.size());
    for (auto const &e : errors) {
      for (auto const &w : e.warnings) {
        OUT() << _FORMAT(
            "%s:%d:%d: # (%s)\n\t\t%s\n",
            fileName.c_str(),
            e.line + 1,
            e.column,
            _FORMAT("%s.%s.%s", e.step.c_str(), e.phase.c_str(), e.key.c_str())
                .c_str(),
            w.c_str());
      }
    }
    if (exit_on_warnings) std::exit(1);
  }

  return algorithms;
}
