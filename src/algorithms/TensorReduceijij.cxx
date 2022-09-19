#include <algorithms/TensorReduceijij.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Sisi4s.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorReduceijij);

void TensorReduceijij::run() {

  checkArgumentsOrDie({"Data", "Out"});

  auto T(getTensorArgument<double>("Data"));
  const int Ni(T->lens[0]), Nj(T->lens[1]);
  const std::vector<int> syms({NS, NS}), lens({Ni, Nj});

  auto t(new Tensor<double>(2, lens.data(), syms.data(), *Sisi4s::world));

  (*t)["ij"] = (*T)["ijij"];

  allocatedTensorArgument<double>("Out", t);

  LOG(1, "TensorReduceijij") << "Ni: "<< Ni << std::endl;
  LOG(1, "TensorReduceijij") << "Nj: "<< Nj << std::endl;

}
