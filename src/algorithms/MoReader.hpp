#ifndef _NWCHEM_MOVECS_READER_DEFINED
#define _NWCHEM_MOVECS_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <vector>
#include <map>
#include <string>

namespace sisi4s {

DEFINE_ALGORITHM_HEADER(

    MoReader,

    static std::map<std::string, std::map<std::string, std::string>>
        DEFAULT_SCALINGS,
    DEFAULT_REORDER;

    static std::vector<std::string> BACKENDS;
    const std::string PSI4 = "psi4",
    NWCHEM = "nwchem",
    TURBOMOLE = "turbomole";);

} // namespace sisi4s

#endif
