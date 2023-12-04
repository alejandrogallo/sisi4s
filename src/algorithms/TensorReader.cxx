#include <fstream>
#include <vector>
#include <string>

#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>
#include <algorithms/TensorReader.hpp>

namespace sisi4s {

template <typename F>
static Tensor<F> *
read(Algorithm &alg, const std::string mode, const std::string &name) {
  Tensor<F> *A;
  if (mode == "binary") {
    std::string fileName(alg.getTextArgument("file", name + ".bin"));
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    A = TensorIo::readBinary<F>(fileName);
  } else {
    std::string fileName(alg.getTextArgument("file", name + ".dat").c_str());
    std::string delimiter(alg.getTextArgument("delimiter", " "));
    int64_t bufferSize(
        alg.getIntegerArgument("bufferSize", 128l * 1024 * 1024));
    A = TensorIo::readText<F>(fileName, delimiter, bufferSize);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  }
  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << name;

  int64_t indexCount(1);
  for (int dim(0); dim < A->order; ++dim) { indexCount *= A->lens[dim]; }
  EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;

  return A;
}

IMPLEMENT_ALGORITHM(TensorReader) {
  std::string name(getArgumentData("Data")->getName());

  // make sure all processes start reading the file at the same time in case
  // it has been modified before
  MPI_Barrier(Sisi4s::world->comm);

  const int64_t precision(getIntegerArgument("precision", 64));
  const std::string mode = getTextArgument("mode", "text");
  if (getIntegerArgument("complexVersion", 0) == 1) {
    LOG(0, "TensorReader") << "Reading complex tensor" << std::endl;
    switch (precision) {
    case 64:
      allocatedTensorArgument<sisi4s::complex>(
          "Data",
          read<sisi4s::complex>(*this, mode, name));
      break;
    }
  } else {
    switch (precision) {
    case 64:
      allocatedTensorArgument<Float64>("Data",
                                       read<Float64>(*this, mode, name));
      break;
    }
  }
}

} // namespace sisi4s
