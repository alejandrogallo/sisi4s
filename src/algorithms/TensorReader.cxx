#include <fstream>
#include <vector>
#include <string>

#include <util/TensorIo.hpp>
#include <Step.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {

DEFSPEC(TensorReader,
        SPEC_IN({"file", SPEC_VALUE_DEF("Input file", std::string, "")},
                {"delimiter", SPEC_VALUE_DEF("Delimiter", std::string, " ")},
                {"precision", SPEC_VALUE_DEF("precision", int64_t, 64)},
                {"bufferSize",
                 SPEC_VALUE_DEF("bufferSize", int64_t, 128l * 1024 * 1024)},
                {"complex",
                 SPEC_VALUE_DEF("If reading in a complex tensor", bool, false)},
                {"mode",
                 SPEC_ONE_OF("The mode for a reader",
                             std::string,
                             "binary",
                             "text")}),
        SPEC_OUT({"Data", SPEC_VAROUT("Data out", CTF::Tensor<double> *)}));

template <typename F>
static Tensor<F> *
read(Arguments &in, const std::string mode, const std::string &name) {
  Tensor<F> *A;
  if (mode == "binary") {
    std::string fileName(in.get<std::string>("file"));
    fileName = fileName.size() ? fileName : name + ".bin";
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    A = TensorIo::readBinary<F>(fileName);
  } else {
    std::string fileName(in.get<std::string>("file"));
    fileName = fileName.size() ? fileName : name + ".dat";
    std::string delimiter(in.get<std::string>("delimiter"));
    int64_t bufferSize(in.get<int64_t>("bufferSize"));
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

DEFSTEP(TensorReader) {

  std::string name = out.get_var("Data");

  // make sure all processes start reading the file at the same time in case
  // it has been modified before
  MPI_Barrier(Sisi4s::world->comm);

  const int64_t precision(in.get<int64_t>("precision"));
  const std::string mode = in.get<std::string>("mode");
  if (in.get<bool>("complex")) {
    LOG(0, "TensorReader") << "Reading complex tensor" << std::endl;
    switch (precision) {
    case 64:
      out.set<CTF::Tensor<sisi4s::complex> *>(
          "Data",
          read<sisi4s::complex>(in, mode, name));
      break;
    }
  } else {
    switch (precision) {
    case 64:
      out.set<CTF::Tensor<Float64> *>("Data", read<Float64>(in, mode, name));
      break;
    }
  }
}

} // namespace sisi4s
