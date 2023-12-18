#include <Step.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

/**
 * \brief Testing environement
 */

DEFSPEC(TensorContraction,
        SPEC_IN({"alpha", SPEC_VALUE_DEF("TODO: DOC", double, 1.0)},
                {"beta", SPEC_VALUE_DEF("TODO: DOC", double, 0.0)},
                {"AIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"BIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"ResultIndex", SPEC_VALUE("TODO: DOC", std::string)},
                {"A", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
                {"B", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
                {"Result", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT());

DEFSTEP(TensorContraction) {
  Tensor<double> *A(in.get<Tensor<double> *>("A"));
  Tensor<double> *B(in.get<Tensor<double> *>("B"));
  Tensor<double> *C(in.get<Tensor<double> *>("Result"));
  C->contract(in.get<double>("alpha", 1.0),
              *A,
              in.get<std::string>("AIndex").c_str(),
              *B,
              in.get<std::string>("BIndex").c_str(),
              in.get<double>("beta", 0.0),
              in.get<std::string>("ResultIndex").c_str());
}
