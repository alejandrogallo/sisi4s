#include <algorithms/Delete.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(Delete, SPEC_IN(), SPEC_OUT());

IMPLEMENT_ALGORITHM(Delete) {
  Data *data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    delete data;
    // remention it in case it will be written to in the future
    new Data(dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}

void Delete::dryRun() {
  Data *data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    delete data;
    // remention it in case it will be written to in the future
    new Data(dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}
