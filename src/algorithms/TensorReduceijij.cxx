#include <algorithms/TensorReduceijij.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Sisi4s.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorReduceijij) {}

DEFSPEC(TensorReduceijij,
        SPEC_IN({"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"Out", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(TensorReduceijij) {

  auto T(in.get<Tensor<double> *>("Data"));
  const int Ni(T->lens[0]), Nj(T->lens[1]);
  const std::vector<int> syms({NS, NS}), lens({Ni, Nj});

  auto t(new Tensor<double>(2, lens.data(), syms.data(), *Sisi4s::world));

  (*t)["ij"] = (*T)["ijij"];

  out.set<Tensor<double> *>("Out", t);

  LOG(1, "TensorReduceijij") << "Ni: " << Ni << std::endl;
  LOG(1, "TensorReduceijij") << "Nj: " << Nj << std::endl;
}
