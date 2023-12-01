#include <algorithms/TensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <iomanip>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

namespace sisi4s {

static const std::type_info &check_type(Data *tensor_data) {
  TensorData<double> *real_tensor_data(
      dynamic_cast<TensorData<double> *>(tensor_data));
  if (real_tensor_data) return typeid(double);
  TensorData<sisi4s::complex> *imag_tensor_data(
      dynamic_cast<TensorData<sisi4s::complex> *>(tensor_data));
  if (imag_tensor_data) return typeid(sisi4s::complex);
  throw "Could not detect type of integrals in TensorAntisymmetrizer";
}

IMPLEMENT_ALGORITHM(TensorWriter) {
  std::string dataName(getArgumentData("Data")->getName());
  // do some switch case if in the future you want to implement
  // precission
  if (typeid(double) == check_type(getArgumentData("Data"))) {
    LOG(1, "TensorWriter") << "Writing real tensor" << std::endl;
    write<Float64>(dataName);
  } else {
    LOG(1, "TensorWriter") << "Writing complex tensor" << std::endl;
    write<sisi4s::complex>(dataName);
  }
}

template <typename F>
void TensorWriter::write(const std::string &name) {
  Tensor<F> *A(getTensorArgument<F>("Data"));
  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << name;
  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    // write binary
    std::string fileName(getTextArgument("file", name + ".bin"));
    TensorIo::writeBinary<F>(fileName, *A);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  } else {
    // write text
    std::string fileName(getTextArgument("file", name + ".dat"));
    std::string rowIndexOrder(getTextArgument("rowIndexOrder", ""));
    std::string columnIndexOrder(getTextArgument("columnIndexOrder", ""));
    std::string delimiter(getTextArgument("delimiter", " "));
    TensorIo::writeText<F>(fileName,
                           *A,
                           rowIndexOrder,
                           columnIndexOrder,
                           delimiter);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  }

  int64_t indexCount(1);
  for (int dim(0); dim < A->order; ++dim) { indexCount *= A->lens[dim]; }
  EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;
}

} // namespace sisi4s
