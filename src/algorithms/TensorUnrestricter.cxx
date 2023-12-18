#include <algorithms/TensorUnrestricter.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <util/Tensor.hpp>
#include <numeric>
#include <map>

using namespace sisi4s;

template <typename F>
static Tensor<F> *unrestrictTensor(Tensor<F> &tensor) {

  // check order of tensor
  const std::set<int> _supportedLens({1, 2, 4});
  if (_supportedLens.count(tensor.order) == 0) {
    throw new EXCEPTION("Tensor size not supported");
  }

  //      lens
  // {N0, N1, ... , Nn} for a n-dimensional tensor
  std::vector<int> lens(tensor.lens, tensor.lens + tensor.order);
  std::vector<int> syms(tensor.sym, tensor.sym + tensor.order);
  std::for_each(lens.begin(), lens.end(), [](int &i) { i *= 2; });

  auto result(new Tensor<F>(tensor.order,
                            lens.data(),
                            syms.data(),
                            *Sisi4s::world,
                            ("u" + std::string(tensor.name)).c_str()));

  // vector of vectors of int
  // { MapEquivalenceForLens0, MapEquivalenceForLens1, ... }
  std::vector<std::vector<int>> permuter(tensor.order);
  std::vector<int *> permuterIds(tensor.order);

  // resize the equivalence mappings
  for (unsigned int i(0); i < lens.size(); i++) {
    permuter[i].resize(tensor.lens[i]);
  }

  std::vector<std::function<void(int &)>> transformers({
      // gerade
      [](int &i) { i = 2 * i; },
      // ungerade
      [](int &i) { i = 2 * i + 1; },
  });

  if (tensor.order == 1 || tensor.order == 2)
    for (auto &transformer : transformers) {
      for (unsigned int i(0); i < lens.size(); i++) {
        auto &mapping = permuter[i];
        std::iota(mapping.begin(), mapping.end(), 0);
        std::for_each(mapping.begin(), mapping.end(), transformer);
        permuterIds[i] = mapping.data();
      }
      result->permute(1.0, tensor, permuterIds.data(), 1.0);
    }

  if (tensor.order == 4)
    for (auto &transLeft : transformers) {
      for (auto &transRight : transformers) {
        for (auto &i : std::vector<int>({0, 2})) {
          for (auto &j : std::vector<int>({1, 3})) {
            auto &mappingLeft = permuter[i];
            auto &mappingRight = permuter[j];

            std::iota(mappingLeft.begin(), mappingLeft.end(), 0);
            std::iota(mappingRight.begin(), mappingRight.end(), 0);

            std::for_each(mappingLeft.begin(), mappingLeft.end(), transLeft);
            std::for_each(mappingRight.begin(), mappingRight.end(), transRight);

            permuterIds[i] = mappingLeft.data();
            permuterIds[j] = mappingRight.data();
          }
        }
        result->permute(1.0, tensor, permuterIds.data(), 1.0);
      }
    }

  return result;
}

IMPLEMENT_EMPTY_DRYRUN(TensorUnrestricter) {}

DEFSPEC(TensorUnrestricter,
        SPEC_IN({"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()}),
        SPEC_OUT({"Out", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(TensorUnrestricter) {

  if (in.is_of_type<Tensor<double> *>("Data")) {
    out.set<Tensor<double> *>(
        "Out",
        unrestrictTensor<double>(*in.get<Tensor<double> *>("Data")));
  } else {
    out.set<Tensor<sisi4s::complex> *>(
        "Out",
        unrestrictTensor<sisi4s::complex>(
            *in.get<Tensor<sisi4s::complex> *>("Data")));
  }
}
