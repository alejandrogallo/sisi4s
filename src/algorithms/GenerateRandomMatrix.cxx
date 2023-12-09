#include <algorithms/GenerateRandomMatrix.hpp>
#include <math/Complex.hpp>
#include <math/RandomTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(
    GenerateRandomMatrix,
    SPEC_IN({"n", SPEC_VALUE_DEF("Number of rows", int64_t, 5)},
            {"m", SPEC_VALUE_DEF("Number of columns", int64_t, 5)},
            {"complex", SPEC_VALUE_DEF("Create a complex tensor", bool, false)},
            {"symmetry",
             SPEC_ONE_OF("Symmetry of the matrix",
                         std::string,
                         "none",
                         "hermitian",
                         "anti",
                         "hollow")}),
    SPEC_OUT({"Result", SPEC_VAROUT("The CTF tensor name", Tensor<double> *)}));

template <typename F>
static void generate(std::vector<int> const &lens,
                     std::vector<int> const &syms,
                     Arguments &out) {
  Tensor<F> *C =
      new Tensor<F>(2, lens.data(), syms.data(), *Sisi4s::world, "C");
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  out.set_force<Tensor<F> *>("Result", C);
}

IMPLEMENT_ALGORITHM(GenerateRandomMatrix) {
  int sym = NS;
  std::string symmetry(in.get<std::string>("symmetry"));
  if (symmetry == "hermitian") {
    sym = SY;
  } else if (symmetry == "anti") {
    sym = AS;
  } else if (symmetry == "hollow") {
    sym = SH;
  }
  std::vector<int> lens = {in.get<int64_t>("n"), in.get<int64_t>("m")},
                   syms(2, NS);
  std::cout << lens[0] << std::endl;
  in.get<bool>("complex") ? generate<sisi4s::complex>(lens, syms, out)
                          : generate<double>(lens, syms, out);
}

IMPLEMENT_EMPTY_DRYRUN(GenerateRandomMatrix) {}
