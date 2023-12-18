#include <Step.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <sstream>
#include <Sisi4s.hpp>
#include <util/Parsing.hpp>

using namespace sisi4s;

template <typename F>
static std::string _showVector(const std::vector<F> &v) {
  std::stringstream s;
  for (const auto &i : v) s << i << ", ";
  return "{" + s.str() + "}";
}

DEFSPEC(TensorSlicer,
        SPEC_IN({"Begin", SPEC_VALUE("TODO: DOC", std::vector<int>)},
                {"End", SPEC_VALUE("TODO: DOC", std::vector<int>)},
                {"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"Out", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

DEFSTEP(TensorSlicer) {

  const std::vector<int> begin(in.get<std::vector<int>>("Begin")),
      end(in.get<std::vector<int>>("End"));

  LOG(0, "TensorSlicer") << "begin: " << _showVector(begin) << std::endl;
  LOG(0, "TensorSlicer") << "end: " << _showVector(end) << std::endl;

  auto T(in.get<Tensor<double> *>("Data"));

  LOG(0, "TensorSlicer") << "slicing..." << std::endl;
  auto t(new Tensor<double>(T->slice(begin.data(), end.data())));

  LOG(0, "TensorSlicer") << "lens: "
                         << _showVector(
                                std::vector<int>(t->lens, t->lens + t->order))
                         << std::endl;
  LOG(0, "TensorSlicer") << "sym: "
                         << _showVector(
                                std::vector<int>(t->sym, t->sym + t->order))
                         << std::endl;

  out.set<Tensor<double> *>("Out", t);
}
