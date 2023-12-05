#include <algorithms/TensorGetMax.hpp>
#include <util/Log.hpp>
#include <iostream>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorGetMax);

IMPLEMENT_EMPTY_DRYRUN(TensorGetMax) {}

void TensorGetMax::run() {
  auto tensor(getTensorArgument<double>("Data"));
  double max{tensor->norm_infty()};
  LOG(1, "TensorGetMax") << tensor->get_name() << ":"
                         << "max: " << max << std::endl;
}
