#include <algorithms/TensorUnrestricter.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <map>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorUnrestricter);


CTF::Tensor<double>*
unrestrictTensor(CTF::Tensor<double> &tensor) {

  std::vector<int> lens(tensor.lens, tensor.lens + tensor.order);
  std::vector<int> syms(tensor.sym, tensor.sym + tensor.order);
  std::for_each(lens.begin(), lens.end(), [](int& i) {i *= 2;});

  auto result(new CTF::Tensor<double>(
    tensor.order, lens.data(), syms.data(), *Cc4s::world,
    ("u" + std::string(tensor.name)).c_str() ));
  std::vector< std::vector<int> > permuter(tensor.order);
  std::vector< int* > permuterIds(tensor.order);

  std::vector< std::function<void(int&)> > transformers({
      [](int &i) { i = 2*i; },
      [](int &i) { i = 2*i + 1; },
  });

  // UpUp part
  for (auto& transformer: transformers) {
    for (unsigned int i(0) ; i < lens.size() ; i++) {
      auto& mapping = permuter[i];
      mapping.resize(tensor.lens[i]);
      std::iota(mapping.begin(), mapping.end(), 0);
      std::for_each(mapping.begin(), mapping.end(), transformer);
      permuterIds[i] = mapping.data();
    }
    result->permute(1.0, tensor, permuterIds.data(), 1.0);
  }

  return result;
}

void TensorUnrestricter::run() {
  allocatedTensorArgument<double>(
    "Out",
    unrestrictTensor(
      *getTensorArgument<double>("Data")));
}
