#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/LineNumberStream.hpp>
#include <string>
#include <vector>

namespace sisi4s {

enum InputFileFormat { YAML, CC4S };

/**
 * \brief Parser for sisi4s files specifying the calculation plan, i.e.
 * which algorithms to use in which order.
 */
template <InputFileFormat fmt>
class InputFileParser {

  /**
   * \brief Creates a new interpreter for a sisi4s file of the given name.
   * Upon creation the file will be openend but not yet read.
   */
  InputFileParser(std::string const &fileName);

  /**
   * \brief Parses the sisi4s algorithms contained in the stream.
   * This method must be called with the same stream content on all processes.
   */
  virtual std::vector<Algorithm *> parse() = 0;
};

template <>
class InputFileParser<InputFileFormat::YAML> {
public:
  InputFileParser(std::string const &fileName, bool exit_on_warnings);
  ~InputFileParser();
  std::vector<Algorithm *> parse();

  bool exit_on_warnings = false;

protected:
  std::string fileName;
};

} // namespace sisi4s

#endif
