#include <algorithms/ComplexTensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <iomanip>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(ComplexTensorWriter) {}


DEFSPEC(ComplexTensorWriter,
        SPEC_IN({"columnIndexOrder",
                 SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
                {"delimiter", SPEC_VALUE_DEF("TODO: DOC", std::string, " ")},
                {"file",
                 SPEC_VALUE_DEF("TODO: DOC", std::string, dataName + ".bin")},
                {"file",
                 SPEC_VALUE_DEF("TODO: DOC", std::string, dataName + ".dat")},
                {"mode", SPEC_VALUE_DEF("TODO: DOC", std::string, "text")},
                {"rowIndexOrder", SPEC_VALUE_DEF("TODO: DOC", std::string, "")},
                {"Data", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ComplexTensorWriter) {
  Tensor<complex> *A(in.get<Tensor<complex> *>("Data"));
  std::string dataName(getArgumentData("Data")->getName());
  A->set_name(dataName.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << dataName;
  std::string mode(in.get<std::string>("mode", "text"));
  if (mode == "binary") {
    // write binary
    std::string fileName(in.get<std::string>("file", dataName + ".bin"));
    TensorIo::writeBinary<complex>(fileName, *A);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  } else {
    // write text
    std::string fileName(in.get<std::string>("file", dataName + ".dat"));
    std::string rowIndexOrder(in.get<std::string>("rowIndexOrder", ""));
    std::string columnIndexOrder(in.get<std::string>("columnIndexOrder", ""));
    std::string delimiter(in.get<std::string>("delimiter", " "));
    TensorIo::writeText<complex>(fileName,
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
