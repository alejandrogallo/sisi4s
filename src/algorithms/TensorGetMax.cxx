#include <Step.hpp>

using namespace sisi4s;

DEFSPEC(TensorGetMax,
        SPEC_IN({"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT());

DEFSTEP(TensorGetMax) {
  auto tensor(in.get<Tensor<double> *>("Data"));
  double max{tensor->norm_infty()};
  LOG(1, "TensorGetMax") << tensor->get_name() << ":"
                         << "max: " << max << std::endl;
}
