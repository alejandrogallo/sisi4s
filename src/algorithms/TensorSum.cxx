#include <algorithms/TensorSum.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorSum) {}


DEFSPEC(TensorSum,
        SPEC_IN({"AFactor", SPEC_VALUE("TODO: DOC", double)},
                {"BFactor", SPEC_VALUE("TODO: DOC", double)},
                {"AIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"BIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"ResultIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"A", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
                {"B", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
                {"Result", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(TensorSum) {
  Tensor<double> *A(in.get<Tensor<double> *>("A"));
  Tensor<double> *B(in.get<Tensor<double> *>("B"));
  Tensor<double> *C(in.get<Tensor<double> *>("Result"));
  (*C)[in.get<std::string>("ResultIndex").c_str()] =
      in.get<double>("AFactor") * (*A)[in.get<std::string>("AIndex").c_str()]
      + in.get<double>("BFactor") * (*B)[in.get<std::string>("BIndex").c_str()];
}
