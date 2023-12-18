#include <fstream>
#include <iomanip>

#include <Step.hpp>
#include <util/TensorIo.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {

DEFSPEC(
    TensorWriter,
    SPEC_IN(
        {"Data", SPEC_VARIN("Input tensor", CTF::Tensor<double> *)},
        {"file", SPEC_VALUE_DEF("File to write it to", std::string, "")},
        {"rowIndexOrder", SPEC_VALUE_DEF("Row index order", std::string, " ")},
        {"columnIndexOrder",
         SPEC_VALUE_DEF("Column index order", std::string, " ")},
        {"delimiter", SPEC_VALUE_DEF("Delimiter", std::string, " ")},
        {"precision", SPEC_VALUE_DEF("precision", int64_t, 64)},
        {"mode",
         SPEC_ONE_OF("The mode for a reader", std::string, "binary", "text")}),
    SPEC_OUT({"Data", SPEC_VAROUT("Data out", CTF::Tensor<double> *)}));

DEFSTEP(TensorWriter) {

  const bool binary_p = in.get<std::string>("mode") == "binary";

  const std::string /**/

      dataName = in.get_var("Data"),
      def_file = in.get<std::string>("file"),
      fileName =
          def_file.size() ? def_file : dataName + (binary_p ? ".bin" : ".dat"),

      rowIndexOrder(in.get<std::string>("rowIndexOrder")),
      columnIndexOrder(in.get<std::string>("columnIndexOrder")),
      delimiter(in.get<std::string>("delimiter"));

  if (in.is_of_type<Tensor<double> *>("Data")) {
    LOG(1, "TensorWriter") << "Writing real tensor to " << fileName
                           << std::endl;
    TensorIo::do_write<double>(dataName,
                               fileName,
                               in.get<Tensor<double> *>("Data"),
                               binary_p,
                               rowIndexOrder,
                               columnIndexOrder,
                               delimiter);
  } else if (in.is_of_type<Tensor<sisi4s::complex> *>("Data")) {
    LOG(1, "TensorWriter") << "Writing complex tensor to " << fileName
                           << std::endl;
    TensorIo::do_write<sisi4s::complex>(
        dataName,
        fileName,
        in.get<Tensor<sisi4s::complex> *>("Data"),
        binary_p,
        rowIndexOrder,
        columnIndexOrder,
        delimiter);
  } else {
    LOG(1, "TensorWriter")
        << "ERROR: I do not know how to write out the tensor " << dataName
        << std::endl;
  }
}

} // namespace sisi4s
