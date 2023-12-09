#include <algorithms/Exit.hpp>
#include <cstdlib>

using namespace sisi4s;

DEFSPEC(Exit, SPEC_IN(), SPEC_OUT());

IMPLEMENT_ALGORITHM(Exit) { std::exit(0); }

IMPLEMENT_EMPTY_DRYRUN(Exit) {}
