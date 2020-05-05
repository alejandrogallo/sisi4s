#ifndef _NWCHEM_MOVECS_READER_DEFINED
#define _NWCHEM_MOVECS_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <vector>

namespace cc4s {

  class NwchemMovecsReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(NwchemMovecsReader);
    NwchemMovecsReader(
      std::vector<Argument> const &argumentList
    ): Algorithm(argumentList) {}
    ~NwchemMovecsReader() {}
    virtual void run();
  };

}

#endif

