#include <algorithms/TensorSlicer.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <sstream>
#include <Sisi4s.hpp>
#include <util/Parsing.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(TensorSlicer) {}

template <typename F>
std::string _showVector(const std::vector<F> &v) {
  std::stringstream s;
  for (const auto &i : v) s << i << ", ";
  return "{" + s.str() + "}";
}

DEFSPEC(TensorSlicer,
        SPEC_IN({"Begin", SPEC_VALUE("TODO: DOC", std::string)},
                {"End", SPEC_VALUE("TODO: DOC", std::string)},
                {"Data", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"Out", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(TensorSlicer) {

  std::vector<int> begin(pars::parseVector<int>(in.get<std::string>("Begin"))),
      end(pars::parseVector<int>(in.get<std::string>("End")));

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
