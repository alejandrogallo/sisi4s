#include<algorithms/Exit.hpp>
#include<cstdlib>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(Exit);

void
cc4s::Exit::run() {
    std::exit(0);
}
