#include<algorithms/Exit.hpp>
#include<cstdlib>

using namespace sisi4s;
ALGORITHM_REGISTRAR_DEFINITION(Exit);

void
sisi4s::Exit::run() {
    std::exit(0);
}
