#include <algorithms/TensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <iomanip>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

namespace sisi4s {

const std::type_info &TensorWriter::check_type(Data *tensor_data) {
  TensorData<double> *real_tensor_data(
      dynamic_cast<TensorData<double> *>(tensor_data));
  if (real_tensor_data) return typeid(double);
  TensorData<sisi4s::complex> *imag_tensor_data(
      dynamic_cast<TensorData<sisi4s::complex> *>(tensor_data));
  if (imag_tensor_data) return typeid(sisi4s::complex);
  throw "Could not detect type of integrals in TensorAntisymmetrizer";
}

IMPLEMENT_ALGORITHM(TensorWriter) {
  const bool binary_p = getTextArgument("mode", "text") == "text";

  const std::string /**/

      dataName = getArgumentData("Data")->getName(),
      fileName = dataName + (binary_p ? ".bin" : ".dat"),

      rowIndexOrder(getTextArgument("rowIndexOrder", "")),
      columnIndexOrder(getTextArgument("columnIndexOrder", "")),
      delimiter(getTextArgument("delimiter", " "));

  if (typeid(double) == TensorWriter::check_type(getArgumentData("Data"))) {
    LOG(1, "TensorWriter") << "Writing real tensor" << std::endl;
    TensorWriter::write<Float64>(dataName,
                                 fileName,
                                 getTensorArgument<Float64>("Data"),
                                 binary_p,
                                 rowIndexOrder,
                                 columnIndexOrder,
                                 delimiter);
  } else {
    LOG(1, "TensorWriter") << "Writing complex tensor" << std::endl;
    TensorWriter::write<sisi4s::complex>(
        dataName,
        fileName,
        getTensorArgument<sisi4s::complex>("Data"),
        binary_p,
        rowIndexOrder,
        columnIndexOrder,
        delimiter);
  }
}

template <typename F>
void TensorWriter::write(const std::string &name,
                         const std::string fileName,
                         Tensor<F> *A,
                         const bool binary_p,
                         const std::string rowIndexOrder,
                         const std::string columnIndexOrder,
                         const std::string delimiter) {

  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << name;
  if (binary_p) {
    TensorIo::writeBinary<F>(fileName, *A);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  } else {
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
