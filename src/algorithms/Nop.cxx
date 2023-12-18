#include <Step.hpp>
using namespace sisi4s;

DEFSPEC(Nop, SPEC_IN(), SPEC_OUT());
DEFSTEP(Nop) {}

DEFSPEC(Lala,
        SPEC_IN({"name", SPEC_VALUE("name", std::string)->require()}),
        SPEC_OUT());
DEFSTEP(Lala) {
  std::cout << "Lalal : " << in.get<std::string>("name") << std::endl;
}
