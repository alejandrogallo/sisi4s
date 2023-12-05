#ifndef FCIDUMP_READER_DEFINED
#define FCIDUMP_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    FcidumpReader,

    struct FcidumpHeader {
      size_t norb;
      size_t nelec;
      size_t uhf;
      size_t ms2;
    };

    FcidumpHeader parseHeader(const std::string &););
} // namespace sisi4s

#endif
