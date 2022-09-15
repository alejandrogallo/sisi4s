#ifndef _ONE_BODYROTATIONDEFINED
#define _ONE_BODYROTATIONDEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

  class OneBodyRotation: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(OneBodyRotation);
    OneBodyRotation(
      std::vector<Argument> const &argumentList
    ): Algorithm(argumentList) {}
    ~OneBodyRotation() {}
    virtual void run();
  };
}

#endif

