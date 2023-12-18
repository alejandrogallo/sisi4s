#include <cstdlib>

#include <Step.hpp>

using namespace sisi4s;

DEFSPEC(Exit, SPEC_IN(), SPEC_OUT());
DEFSTEP(Exit) { std::exit(0); }
