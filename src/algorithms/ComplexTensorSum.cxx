#include <algorithms/ComplexTensorSum.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

/**
 * \brief Testing environement
 */

DEFSPEC(ComplexTensorSum,
        SPEC_IN({"AFactor", SPEC_VALUE("TODO: DOC", double)},
                {"BFactor", SPEC_VALUE("TODO: DOC", double)},
                {"AIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"BIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"ResultIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"A", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"B", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"Result", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ComplexTensorSum) {

  Tensor<complex> *A(in.get<Tensor<complex> *>("A"));
  Tensor<complex> *B(in.get<Tensor<complex> *>("B"));
  Tensor<complex> *C(in.get<Tensor<complex> *>("Result"));
  (*C)[in.get<std::string>("ResultIndex").c_str()] =
      in.get<double>("AFactor") * (*A)[in.get<std::string>("AIndex").c_str()]
      + in.get<double>("BFactor") * (*B)[in.get<std::string>("BIndex").c_str()];
}
