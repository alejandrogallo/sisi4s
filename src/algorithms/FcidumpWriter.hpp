#ifndef FCIDUMP_WRITER_DEFINED
#define FCIDUMP_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

class FcidumpWriter : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(FcidumpWriter);
  FcidumpWriter(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList){};
  ~FcidumpWriter(){};
  virtual void run();
};
} // namespace sisi4s

#endif
