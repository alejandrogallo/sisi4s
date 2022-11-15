/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <Data.hpp>
#include <string>
#include <vector>
#include <sstream>
#include <util/Tensor.hpp>

namespace sisi4s {
  class Argument {
  public:
    Argument(
      std::string const &name_
    ): name(name_), data(name_) {
    }
    Argument(
      std::string const &name_, std::string const &data_
    ): name(name_), data(data_) {
    }
    std::string const &getName() const { return name; }
    std::string const &getData() const { return data; }
  protected:
    std::string name, data;
  };

  class Algorithm {
  public:
    Algorithm(std::vector<Argument> const &argumentList);
    virtual ~Algorithm();
    virtual std::string getName() = 0;
    virtual void run() = 0;
    virtual void dryRun();

    bool isArgumentGiven(std::string const &argumentName);
    // retrieving input arguments
    std::string getTextArgument(std::string const &argumentName);
    std::string getTextArgument(
      std::string const &argumentName, std::string const &defaultValue
    );
    bool getBooleanArgument(std::string const &name);
    bool getBooleanArgument(
      std::string const &name, bool const &defaultValue
    );
    int64_t getIntegerArgument(std::string const &argumentName);
    int64_t getIntegerArgument(
      std::string const &argumentName, int64_t const defaultValue
    );
    real getRealArgument(std::string const &argumentName);
    real getRealArgument(
      std::string const &argumentName, real const defaultValue
    );
    template < typename F=real, typename T=Tensor<F> >
    T *getTensorArgument(std::string const &argumentName);
    template < typename F=real, typename C=std::vector<F> >
    C *getContainerArgument(std::string const &argumentName);
    template < typename F=real, typename C=std::vector<F> >
    void allocateContainerArgument(
      std::string const &argumentName, C *container
    );

    // Get all arguments given by the user in the input file
    std::vector<std::string> getGivenArgumentNames() const {
      std::vector<std::string> names;
      for (const auto& p: arguments) { names.push_back(p.first); }
      return names;
    }

    // Check that all the arguments given by the user
    // conform to the args that you're expecting in the algorithm
    void checkArgumentsOrDie(const std::vector<std::string> args) const {
      const std::vector<std::string> actualArguments(getGivenArgumentNames());
      std::stringstream s;
      for (const auto& p: actualArguments) {
        if (std::count(args.begin(), args.end(), p) != 1) {
          s << "Error: parameter (" << p << ") unknown";
          throw new EXCEPTION(s.str());
        }
      }
    }

    // typing, allocating and setting output arguments
    /**
     * \brief Specifies the location of an output tensor data.
     * \param[in] argumentName The argument name as specified in the sisi4s file
     * \param[in] tensor The reference of the tensor data allocated by the
     * caller and later freed by the system if not needed any further.
     * \note
     * often explicit instantiation may be necessary, e.g.
     * \code{.cxx}
     * allocatedTensorArgument<complex>(complexTensor);
     * \endcode
     */
    template < typename F=real, typename T=Tensor<F> >
    void allocatedTensorArgument(
      std::string const &argumentName, T *tensor
    );
    void setRealArgument(std::string const &argumentName, real const value);
    void setIntegerArgument(std::string const &argumentName, int const value);

  protected:
    // type promotions:
    real getRealArgumentFromInteger(IntegerData *data);
    real getRealArgumentFromTensor(TensorData<real> *data);
    template < typename F=real, typename T=Tensor<F> >
    T *getTensorArgumentFromReal(RealData *realData);

    Data *getArgumentData(std::string const &argumentName);
    std::map<std::string, std::string> arguments;
  };

  class AlgorithmFactory {
  public:

    using AlgorithmMap
      = std::map<std::string,
                 std::function<Algorithm *(std::vector<Argument> const &)>
                 >;

    /**
     * \brief Creates an algorithm object of the algorithm type specified
     * by the given name. The given arguments are passed to the algorithm
     * constructor.
     * The instantiated algorithm must be registered using the
     * AlgorithmRegistrar class.
     */
    static Algorithm *create(
      std::string const &name, std::vector<Argument> const &arguments
    ) {
      auto iterator(getAlgorithmMap()->find(name));
      return iterator != getAlgorithmMap()->end() ?
        iterator->second(arguments) : nullptr;
    }

    static std::vector<std::string>
    getAlgorithmNames() {
      std::vector<std::string> result;
      for (auto const& kv: *algorithmMap) result.push_back(kv.first);
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
  Algorithm *createAlgorithm(std::vector<Argument> const &arguments) {
    return new AlgorithmType(arguments);
  }

  /**
   * \brief Class to be statically instantiated by an algorithm to register
   * it in the AlgorithmFactory. Registered algorithms can be instantiated
   * from the sisi4s control language.
   */
  template <typename AlgorithmType>
  class AlgorithmRegistrar: protected AlgorithmFactory {
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
#define ALGORITHM_REGISTRAR_DECLARATION(NAME)     \
  virtual std::string getName() { return #NAME; } \
  static AlgorithmRegistrar<NAME> registrar_
  /**
   * \brief Auxiliary macro defining the algorithm registrar for
   * the algorithm type of the given name. This macro is to be
   * used in the algorithm definition within the .cxx file.
   * Note that name is a symbol name not a string.
   */
#define ALGORITHM_REGISTRAR_DEFINITION(NAME)        \
  AlgorithmRegistrar<NAME> NAME::registrar_(#NAME)

#define IMPLEMENT_ALGORITHM(NAME)                 \
  ALGORITHM_REGISTRAR_DEFINITION(NAME);           \
  void NAME::run()

#define DEFINE_ALGORITHM_HEADER(NAME, ...)        \
  class NAME: public Algorithm {                  \
  public:                                         \
  ALGORITHM_REGISTRAR_DECLARATION(NAME);          \
  NAME(std::vector<Argument> const &argumentList) \
    : Algorithm(argumentList) {}                  \
  ~NAME(){}                                       \
  virtual void run();                             \
  __VA_ARGS__                                     \
  }

}

#endif
