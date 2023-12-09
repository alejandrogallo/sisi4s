#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <algorithm>

#include <NewData.hpp>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <util/Format.hpp>
#include <AlgorithmInputSpec.hpp>

namespace sisi4s {

class Arguments {
public:
  using ArgumentName = std::string;
  using Index = std::string;
  std::map<ArgumentName, Index> arguments;

  static std::string normalize_name(std::string const &name) {
    const std::regex symbols("[-_ ]");
    std::string new_name = std::regex_replace(name, symbols, "");
    std::transform(new_name.begin(),
                   new_name.end(),
                   new_name.begin(),
                   [](unsigned char const &c) { return std::tolower(c); });
    return new_name;
  }

  void push(std::string const &name, std::string const &db_index) {
    arguments[normalize_name(name)] = db_index;
  }

  bool present(std::string const &name) {
    return arguments.find(normalize_name(name)) != arguments.end();
  }

  data::StorePair get_data(std::string const &name) {
    const std::string nname = normalize_name(name);
    auto dataIterator(arguments.find(nname));
    if (dataIterator == arguments.end()) {
      std::stringstream sStream;
      sStream << "Missing argument: " << name;
      //    throw new EXCEPTION(std::stringstream() << "Missing argument: " <<
      //    name);
      throw new EXCEPTION(sStream.str());
    }
    std::string db_index = dataIterator->second;
    // Data *data = Data::get(db_index);
    data::Data *data = data::getraw(db_index);
    if (!data) {
      throw new EXCEPTION(_FORMAT("Missing data in db '%s' for key '%s' ",
                                  db_index.c_str(),
                                  name.c_str()));
    }
    return {db_index, data};
  }

  std::string get_var(std::string const &name) {
    return get_data(normalize_name(name)).second->name;
  }

  template <typename F>
  bool is_of_type(std::string const &name) {
    return data::istype<F>(get_data(normalize_name(name)).first);
  }

  template <typename F>
  F get(std::string const &name) {
    const std::string nname = normalize_name(name);
    F value = *getptr<F>(nname);
    EMIT() << YAML::Key << nname << YAML::Value << value;
    return value;
  }

  template <typename F>
  F *getptr(std::string const &name) {
    const std::string nname = normalize_name(name);
    auto const pair = get_data(nname);
    const std::string db_index = pair.first;
    return data::get<F>(db_index);
  }

  // template <typename F>
  // F get(std::string const &name, F const &def) {
  //   if (!present(name)) {
  //     EMIT() << YAML::Key << name << YAML::Value << def;
  //     return def;
  //   } else {
  //     return get<F>(name);
  //   }
  // }

  template <typename F>
  void set(std::string const &name, F const &value) {
    F *ptr = getptr<F>(normalize_name(name));
    *ptr = value;
  }

  template <typename F>
  void set_force(std::string const &name, F const &value) {
    const std::string nname = normalize_name(name);
    auto const pair = get_data(nname);
    const std::string db_index = pair.first;
    set<F>(nname, value);
    data::declare_type<F>(db_index);
  }
};

class Algorithm {
public:
  using ArgumentMap = std::map<std::string, std::string>;
  Algorithm(Arguments const &in, Arguments const &out);
  virtual ~Algorithm();
  virtual std::string getName() = 0;
  virtual void run() = 0;
  virtual void dryRun();
  std::string note;
  bool fallible = false;
  Arguments in, out;
};

class AlgorithmFactory {
public:
  using AlgorithmMap = std::map<
      std::string,
      std::function<Algorithm *(Arguments const &, Arguments const &)>>;

  /**
   * \brief Creates an algorithm object of the algorithm type specified
   * by the given name. The given arguments are passed to the algorithm
   * constructor.
   * The instantiated algorithm must be registered using the
   * AlgorithmRegistrar class.
   */
  static Algorithm *
  create(std::string const &name, Arguments const &in, Arguments const &out) {
    auto iterator(getAlgorithmMap()->find(name));
    return iterator != getAlgorithmMap()->end() ? iterator->second(in, out)
                                                : nullptr;
  }

  static std::vector<std::string> getAlgorithmNames() {
    std::vector<std::string> result;
    for (auto const &kv : *algorithmMap) result.push_back(kv.first);
    return result;
  }

protected:
  static AlgorithmMap *getAlgorithmMap() {
    return algorithmMap ? algorithmMap : (algorithmMap = new AlgorithmMap);
  }
  static AlgorithmMap *algorithmMap;
};

/**
 * \brief template function creating an instance of the given class.
 */
template <typename AlgorithmType>
Algorithm *createAlgorithm(Arguments const &in, Arguments const &out) {
  return new AlgorithmType(in, out);
}

/**
 * \brief Class to be statically instantiated by an algorithm to register
 * it in the AlgorithmFactory. Registered algorithms can be instantiated
 * from the sisi4s control language.
 */
template <typename AlgorithmType>
class AlgorithmRegistrar : protected AlgorithmFactory {
public:
  /**
   * \brief Constructs the registrating instance. The algorithm type
   * must be given as template argument, the algorithm name as
   * method argument.
   */
  AlgorithmRegistrar(std::string const &name) {
    (*getAlgorithmMap())[name] = &createAlgorithm<AlgorithmType>;
  }
};

/**
 * \brief Auxiliary macro declaring the algorithm registrar for
 * the algorithm type of the given name. This macro is to be
 * used in the algorith declaration within the .hpp file.
 * Note that name is a symbol name not a string.
 */
#define ALGORITHM_REGISTRAR_DECLARATION(NAME)                                  \
  virtual std::string getName() { return #NAME; }                              \
  static sisi4s::AlgorithmRegistrar<NAME> registrar_
/**
 * \brief Auxiliary macro defining the algorithm registrar for
 * the algorithm type of the given name. This macro is to be
 * used in the algorithm definition within the .cxx file.
 * Note that name is a symbol name not a string.
 */
#define ALGORITHM_REGISTRAR_DEFINITION(NAME)                                   \
  sisi4s::AlgorithmRegistrar<NAME> NAME::registrar_(#NAME)

#define IMPLEMENT_ALGORITHM(NAME)                                              \
  ALGORITHM_REGISTRAR_DEFINITION(NAME);                                        \
  void NAME::run()

#define IMPLEMENT_EMPTY_DRYRUN(NAME) void NAME::dryRun()

#define DEFINE_ALGORITHM_HEADER(NAME, ...)                                     \
  class NAME : public Algorithm {                                              \
  public:                                                                      \
    ALGORITHM_REGISTRAR_DECLARATION(NAME);                                     \
    using Algorithm::Algorithm;                                                \
    static AlgorithmInputSpec spec;                                            \
    virtual void run();                                                        \
    virtual void dryRun();                                                     \
    __VA_ARGS__                                                                \
  }

} // namespace sisi4s

#endif
