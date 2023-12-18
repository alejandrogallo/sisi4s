#include <Step.hpp>
#include <math/RandomTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(GenerateRandomTensor,
        SPEC_IN({"dimensions",
                 SPEC_VALUE_DEF("Dimensions of a tensor",
                                std::vector<int>,
                                {5, 5, 5, 5})},
                {"complex",
                 SPEC_VALUE_DEF("Create a complex tensor", bool, false)}, ),
        SPEC_OUT({"Result",
                  SPEC_VAROUT("The CTF tensor name", Tensor<double> *)}));

template <typename F>
static void generate(std::vector<int> const &lens,
                     std::vector<int> const &syms,
                     Arguments &out) {
  Tensor<F> *C =
      new Tensor<F>(lens.size(), lens.data(), syms.data(), *Sisi4s::world);
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  out.set_force<Tensor<F> *>("Result", C);
}

DEFSTEP(GenerateRandomTensor) {

  std::vector<int> lens = in.get<std::vector<int>>("dimensions"),
                   syms(lens.size(), NS);
  in.get<bool>("complex") ? generate<sisi4s::complex>(lens, syms, out)
                          : generate<double>(lens, syms, out);
}
