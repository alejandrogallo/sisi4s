#include <algorithms/ComplexTensorContraction.hpp>
#include <math/Complex.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

/**
 * \brief Testing environement
 */

DEFSPEC(ComplexTensorContraction,
        SPEC_IN({"AIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"BIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"ResultIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"A", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"B", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"Result", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ComplexTensorContraction) {

  Tensor<complex> *A(in.get<Tensor<complex> *>("A"));
  Tensor<complex> *B(in.get<Tensor<complex> *>("B"));
  Tensor<complex> *C(in.get<Tensor<complex> *>("Result"));
  (*C)[in.get<std::string>("ResultIndex").c_str()] =
      (*A)[in.get<std::string>("AIndex").c_str()]
      * (*B)[in.get<std::string>("BIndex").c_str()];
}
