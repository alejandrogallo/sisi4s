#ifndef FCIDUMP_READER_DEFINED
#define FCIDUMP_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

class FcidumpReader : public Algorithm {
public:
  struct FcidumpHeader {
    size_t norb;
    size_t nelec;
    size_t uhf;
    size_t ms2;
  };

  FcidumpHeader parseHeader(const std::string &);

  ALGORITHM_REGISTRAR_DECLARATION(FcidumpReader);
  FcidumpReader(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList){};
  ~FcidumpReader(){};
  virtual void run();
};
} // namespace sisi4s

#endif
