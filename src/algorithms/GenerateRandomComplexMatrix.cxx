#include <algorithms/GenerateRandomComplexMatrix.hpp>
#include <math/RandomTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(GenerateRandomComplexMatrix) {}

/**
 * \brief Testing environement
 */

DEFSPEC(GenerateRandomComplexMatrix,
        SPEC_IN({"m", SPEC_VALUE_DEF("TODO: DOC", int64_t, 200)},
                {"n", SPEC_VALUE_DEF("TODO: DOC", int64_t, 200)},
                {"symmetric",
                 SPEC_VALUE_DEF("TODO: DOC", std::string, "none")}),
        SPEC_OUT({"Result", SPEC_VAROUT("TODO: DOC", Tensor<complex> *)}));

IMPLEMENT_ALGORITHM(GenerateRandomComplexMatrix) {
  int m(in.get<int64_t>("m", 200)), n(in.get<int64_t>("n", 200));
  std::string symmetry(in.get<std::string>("symmetric", "none"));
  int sym(NS);
  if (symmetry == "hermitian") {
    throw new EXCEPTION(
        "Hermitian symmetry of complex tensors not yet supported.");
  }
  Matrix<complex> *C(new Matrix<complex>(m, n, sym, *Sisi4s::world, "C"));
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  out.set<Tensor<complex> *>("Result", C);
}
