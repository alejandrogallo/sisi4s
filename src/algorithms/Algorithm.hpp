#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <algorithm>

#include <NewData.hpp>
#include <Arguments.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <util/Format.hpp>
#include <AlgorithmInputSpec.hpp>

namespace sisi4s {

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
  using AlgorithmMaker =
      std::function<Algorithm *(Arguments const &, Arguments const &)>;
  using AlgorithmMap = std::map<std::string, AlgorithmMaker>;

  static std::string normalize_name(std::string const &name);

  static AlgorithmMap::iterator find_by_name(std::string const &name) {
    std::string const &new_name = AlgorithmFactory::normalize_name(name);
    return get_algorithm_map()->find(new_name);
  }

  /**
   * \brief Creates an algorithm object of the algorithm type specified
   * by the given name. The given arguments are passed to the algorithm
   * constructor.
   * The instantiated algorithm must be registered using the
   * AlgorithmRegistrar class.
   */
  static Algorithm *
  create(std::string const &name, Arguments const &in, Arguments const &out) {
    auto iterator = AlgorithmFactory::find_by_name(name);
    return iterator != get_algorithm_map()->end() ? iterator->second(in, out)
                                                  : nullptr;
  }

  static std::vector<std::string> get_algorithm_names() {
    std::vector<std::string> result;
    for (auto const &kv : *algorithmMap) result.push_back(kv.first);
    return result;
  }

protected:
  static AlgorithmMap *get_algorithm_map() {
    return algorithmMap ? algorithmMap : (algorithmMap = new AlgorithmMap);
  }
  static AlgorithmMap *algorithmMap;
};

/**
 * \brief template function creating an instance of the given class.
 */
template <typename AlgorithmType>
Algorithm *make_algorithm(Arguments const &in, Arguments const &out) {
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
    (*get_algorithm_map())[name] = &make_algorithm<AlgorithmType>;
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
    virtual void run();                                                        \
    virtual void dryRun();                                                     \
    __VA_ARGS__                                                                \
  }

#define DEFSTEP(name)                                                          \
  class name;                                                                  \
  static AlgorithmRegistrar<name> name##_registrar(#name);                     \
  class name : public Algorithm {                                              \
    using Algorithm::Algorithm;                                                \
    virtual std::string getName() { return #name; }                            \
    virtual void dryRun() {}                                                   \
    virtual void run();                                                        \
  };                                                                           \
  void name::run()

} // namespace sisi4s

#endif
