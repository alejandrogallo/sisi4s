#include <algorithms/TensorGetMax.hpp>
#include <util/Log.hpp>
#include <iostream>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorGetMax) {}


DEFSPEC(TensorGetMax,
        SPEC_IN({"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(TensorGetMax) {
  auto tensor(in.get<Tensor<double> *>("Data"));
  double max{tensor->norm_infty()};
  LOG(1, "TensorGetMax") << tensor->get_name() << ":"
                         << "max: " << max << std::endl;
}
