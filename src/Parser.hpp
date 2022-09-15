#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <Data.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/LineNumberStream.hpp>
#include <string>
#include <vector>

namespace sisi4s {

  enum InputFileFormat {
    YAML,
    CC4S
  };

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
    InputFileParser(std::string const& fileName);

    /**
     * \brief Parses the sisi4s algorithms contained in the stream.
     * This method must be called with the same stream content on all processes.
     */
    virtual std::vector<Algorithm*> parse() = 0;
  };

  template <>
  class InputFileParser<InputFileFormat::YAML> {
    public:
      InputFileParser(std::string const &fileName);
      ~InputFileParser();
      std::vector<Algorithm *> parse();
    protected:
      std::string fileName;
  };

  template<>
  class InputFileParser<InputFileFormat::CC4S> {
  public:

    InputFileParser(std::string const &fileName);
    ~InputFileParser();

    std::vector<Algorithm *> parse();

  protected:
    Algorithm *parseAlgorithm();
    std::vector<Argument> parseArguments();
    Argument parseArgument();
    Argument parseImplicitlyNamedArgument();
    Argument parseExplicitlyNamedArgument();
    std::string parseData();
    std::string parseSymbolName();
    Data *parseSymbol();
    TextData *parseText();
    NumericData *parseNumber();
    RealData *parseReal(int64_t const sign, int64_t const integerPart);

    void skipIrrelevantCharacters();
    void skipComment();
    void skipWhiteSpaceCharacters();
    void expectCharacter(char const character);

    LineNumberStream stream;
  };
}

#endif

