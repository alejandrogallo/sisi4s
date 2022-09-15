/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <fstream>
#include <vector>
#include <string>

#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/CTF.hpp>
#include <util/Emitter.hpp>
#include <algorithms/TensorReader.hpp>


template<typename F>
using Tensor = CTF::Tensor<F>;

namespace sisi4s {

  ALGORITHM_REGISTRAR_DEFINITION(TensorReader);

  TensorReader::TensorReader(std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

  TensorReader::~TensorReader() {
  }

  void TensorReader::run() {
    std::string name(getArgumentData("Data")->getName());

    // make sure all processes start reading the file at the same time in case
    // it has been modified before
    MPI_Barrier(Sisi4s::world->comm);

    int64_t precision(getIntegerArgument("precision", 64));
    switch (precision) {
    case 64:
      allocatedTensorArgument<Float64>("Data", read<Float64>(name));
      break;
    }
  }

  template <typename F>
  Tensor<F> *TensorReader::read(const std::string &name) {
    Tensor<F> *A;
    std::string mode(getTextArgument("mode", "text"));
    if (mode == "binary") {
      std::string fileName(getTextArgument("file", name + ".bin"));
      EMIT() << YAML::Key << "file" << YAML::Value << fileName;
      A = TensorIo::readBinary<F>(fileName);
    } else {
      std::string fileName(getTextArgument("file", name + ".dat").c_str());
      std::string delimiter(getTextArgument("delimiter", " "));
      int64_t bufferSize(getIntegerArgument("bufferSize", 128l*1024*1024));
      A = TensorIo::readText<F>(fileName, delimiter, bufferSize);
      EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    }
    A->set_name(name.c_str());
    EMIT() << YAML::Key << "Data"  << YAML::Value << name;

    int64_t indexCount(1);
    for (int dim(0); dim < A->order; ++dim) {
      indexCount *= A->lens[dim];
    }
    EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;

    return A;
  }

}  // namespace sisi4s
