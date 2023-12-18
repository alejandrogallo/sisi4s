#ifndef FCIDUMP_HPP_
#define FCIDUMP_HPP_

namespace sisi4s {

struct FcidumpHeader {
  size_t norb;
  size_t nelec;
  size_t uhf;
  size_t ms2;
};

} // namespace sisi4s

#endif
