#ifndef EXIT_ALGORITHM_DEFINED_HERE
#define EXIT_ALGORITHM_DEFINED_HERE

#include <algorithms/Algorithm.hpp>

namespace cc4s {

  class Exit: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Exit);
    Exit(std::vector<Argument> const &argumentList) :Algorithm(argumentList) {}
    virtual ~Exit() {}
    virtual void run();
  };
}



#endif // EXIT_ALGORITHM_DEFINED_HERE
