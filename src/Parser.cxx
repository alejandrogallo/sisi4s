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
  std::string const& fileName_
): fileName(fileName_) {}

InputFileParser<InputFileFormat::YAML>::~InputFileParser() {}

std::vector<Algorithm *>
InputFileParser<InputFileFormat::YAML>::parse() {
  const YAML::Node nodes = YAML::LoadFile(fileName);
  std::vector<Algorithm *> algorithms;

  for (const YAML::Node& node: nodes) {
    std::string name = node["name"].as<std::string>();
    std::vector<Argument> arguments;

    for (auto const& _inout: {"in", "out"}) {
      YAML::Node const inout = node[_inout];
      for (YAML::const_iterator it = inout.begin(); it != inout.end(); ++it) {
        std::string
          key = it->first.as<std::string>(),
          valueName;
        try {
          int value = it->second.as<int>();
          valueName = (new IntegerData(value))->getName();
        } catch (YAML::TypedBadConversion<int> const& c) {

          try {
            double value = it->second.as<double>();
            valueName = (new RealData(value))->getName();
          } catch (YAML::TypedBadConversion<double> const& c) {

            try {
              std::string value = it->second.as<std::string>();
              if (value.substr(0, 1) == "$") {
                const std::string symbolName = value.substr(1);
                Data *data(Data::get(symbolName));
                valueName = data
                          ? data->getName()
                          : (new Data(symbolName))->getName();
              } else {
                valueName = (new TextData(value))->getName();
              }
            } catch (YAML::TypedBadConversion<std::string> const& c) {
              throw c;
            }

          }

        }

        arguments.push_back(Argument(key, valueName));
      }
    }

    Algorithm *algorithm(AlgorithmFactory::create(name, arguments));
    algorithms.push_back(algorithm);
  }
  return algorithms;
}


///////////////////////////////////////////////////////////////////////////////
// Hummel Parser //////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

InputFileParser<InputFileFormat::CC4S>::InputFileParser(
  std::string const &fileName
): stream(new std::ifstream(fileName.c_str()), fileName) {
  std::ifstream *fileStream(dynamic_cast<std::ifstream *>(stream.getStream()));
  if (!fileStream->is_open()) {
    std::stringstream sStream;
    sStream << "Failed to open file " << fileName;
    throw new EXCEPTION(sStream.str());
  }
}

InputFileParser<InputFileFormat::CC4S>::~InputFileParser() {}

std::vector<Algorithm *> InputFileParser<InputFileFormat::CC4S>::parse() {
  std::vector<Algorithm *> algorithms;
  skipIrrelevantCharacters();
  while (stream.peek() > 0) {
    algorithms.push_back(parseAlgorithm());
    skipIrrelevantCharacters();
  }
  return algorithms;
}

Algorithm* InputFileParser<InputFileFormat::CC4S>::parseAlgorithm() {
  int line(stream.getLine()), column(stream.getColumn());
  // an algorithm starts with the name
  std::string algorithmName(parseSymbolName());
  skipIrrelevantCharacters();
  // next comes its input arguments
  std::vector<Argument> arguments(parseArguments());
  skipIrrelevantCharacters();
  // next comes its output arguments
  std::vector<Argument> outputArguments(parseArguments());
  skipIrrelevantCharacters();
  // the algorithm must end with a period
  expectCharacter('.');
  // currently there is no distinction between input and output arguments
  arguments.insert(
    arguments.end(), outputArguments.begin(), outputArguments.end()
  );
  // create and return an instance of the algorithm
  Algorithm *algorithm(AlgorithmFactory::create(algorithmName, arguments));
  if (!algorithm) {
    std::stringstream sStream;
    sStream << "Unknown algorithm " << algorithmName;
    throw new DetailedException(sStream.str(), stream.getSource(), line,column);
  }
  return algorithm;
}

std::vector<Argument> InputFileParser<InputFileFormat::CC4S>::parseArguments() {
  std::vector<Argument> arguments;
  expectCharacter('[');
  skipIrrelevantCharacters();
  while (stream.peek() != ']') {
    arguments.push_back(parseArgument());
    skipIrrelevantCharacters();
  }
  expectCharacter(']');
  return arguments;
}

Argument InputFileParser<InputFileFormat::CC4S>::parseArgument() {
  if (stream.peek() == '(') return parseExplicitlyNamedArgument();
  else return parseImplicitlyNamedArgument();
}

Argument InputFileParser<InputFileFormat::CC4S>::parseImplicitlyNamedArgument() {
  // TODO: store debug info for later reference in case of errors
  std::string argumentName(parseSymbolName());
  Data *data(Data::get(argumentName));
  if (!data) new Data(argumentName);
  return Argument(argumentName, argumentName);
}

Argument InputFileParser<InputFileFormat::CC4S>::parseExplicitlyNamedArgument() {
  // TODO: store debug info for later reference in case of errors
  // first character must be '('
  stream.get();
  skipIrrelevantCharacters();
  std::string argumentName(parseSymbolName());
  skipIrrelevantCharacters();
  std::string dataName(parseData());
  skipIrrelevantCharacters();
  expectCharacter(')');
  return Argument(argumentName, dataName);
}

std::string InputFileParser<InputFileFormat::CC4S>::parseData() {
  char character(stream.peek());
  if (isalpha(character)) {
    return parseSymbol()->getName();
  } else if (isdigit(character) || character == '+' || character == '-') {
    return parseNumber()->getName();
  } else if (character == '"') {
    return parseText()->getName();
  } else {
    throw new DetailedException(
      "Constant or symbol expression expected",
      stream.getSource(), stream.getLine(), stream.getColumn()
    );
  }
}

std::string InputFileParser<InputFileFormat::CC4S>::parseSymbolName() {
  std::stringstream sStream;
  // the first character is expected to be an alphabetic character
  sStream.put(stream.get());
  char c;
  while (isalpha(c = stream.peek()) || isdigit(c)) {
    sStream.put(stream.get());
  }
  return sStream.str();
}

Data *InputFileParser<InputFileFormat::CC4S>::parseSymbol() {
  std::string symbolName(parseSymbolName());
  Data *data(Data::get(symbolName));
  return data ? data : new Data(symbolName);
}

TextData *InputFileParser<InputFileFormat::CC4S>::parseText() {
  std::stringstream sStream;
  // TODO: parse escape sequences
  // the first character is expected to be a double quote '"'
  stream.get();
  char c;
  while ((c = stream.peek()) > 0 && c != '"') {
    sStream.put(stream.get());
  }
  expectCharacter('"');
  return new TextData(sStream.str());
}

NumericData *InputFileParser<InputFileFormat::CC4S>::parseNumber() {
  // the first character can be a sign
  int64_t sign(1);
  switch (stream.peek()) {
  case '-':
    sign = -1;
  case '+':
    stream.get();
  }
  // the next character must be a digit
  if (!isdigit(stream.peek())) {
    throw new DetailedException(
      "Digit expected",
      stream.getSource(), stream.getLine(), stream.getColumn()
    );
  }
  int64_t integer(stream.get() - '0');
  while (isdigit(stream.peek())) {
    integer *= 10;
    integer += stream.get() - '0';
  }
  if (stream.peek() == '.') return parseReal(sign, integer);
  else return new IntegerData(sign * integer);
}

RealData *InputFileParser<InputFileFormat::CC4S>::parseReal(
  const int64_t sign, const int64_t integerPart
) {
  // the first character is expected to be the decimal point
  stream.get();
  int64_t numerator(0), denominator(1);
  while (isdigit(stream.peek())) {
    numerator *= 10; denominator *= 10;
    numerator += stream.get() - '0';
  }
  // TODO: parse scientific notatoin e-1
  return new RealData(
    sign * (integerPart + double(numerator) / denominator)
  );
}

void InputFileParser<InputFileFormat::CC4S>::skipIrrelevantCharacters() {
  skipWhiteSpaceCharacters();
  while (stream.peek() == '%') {
    skipComment();
    skipWhiteSpaceCharacters();
  }
}

void InputFileParser<InputFileFormat::CC4S>::skipComment() {
  char c;
  while ((c = stream.get()) > 0 && c != '\n');
}

void InputFileParser<InputFileFormat::CC4S>::skipWhiteSpaceCharacters() {
  while (isspace(stream.peek())) stream.get();
}

void InputFileParser<InputFileFormat::CC4S>::expectCharacter(char const expectedCharacter) {
  char character(stream.peek());
  if (character != expectedCharacter) {
    std::stringstream sStream;
    sStream <<
      "Expected '" << expectedCharacter << "', got '" << character << "'";
    throw new DetailedException(
      sStream.str(), stream.getSource(), stream.getLine(), stream.getColumn()
    );
  }
  stream.get();
}

