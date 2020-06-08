#ifndef _NWCHEM_MOVECS_READER_DEFINED
#define _NWCHEM_MOVECS_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <map>
#include <string>

namespace cc4s {

  class NwchemMovecsReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(NwchemMovecsReader);
    NwchemMovecsReader(
      std::vector<Argument> const &argumentList
    ): Algorithm(argumentList) {}
    ~NwchemMovecsReader() {}
    virtual void run();

    static std::map<std::string, std::map<std::string, std::string>>
      DEFAULT_SCALINGS, DEFAULT_REORDER;

    static std::vector<std::string> BACKENDS;

  };

}

#endif

