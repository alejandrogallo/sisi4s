/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream>
#include <iomanip>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>

namespace sisi4s {

  IMPLEMENT_ALGORITHM(TensorWriter) {
    std::string dataName(getArgumentData("Data")->getName());
    // do some switch case if in the future you want to implement
    // precission
    write<Float64>(dataName);
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
      TensorIo::writeText<F>(
        fileName, *A, rowIndexOrder, columnIndexOrder, delimiter
      );
      EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    }

    int64_t indexCount(1);
    for (int dim(0); dim < A->order; ++dim) {
      indexCount *= A->lens[dim];
    }
    EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;

  }

}  // namespace sisi4s
