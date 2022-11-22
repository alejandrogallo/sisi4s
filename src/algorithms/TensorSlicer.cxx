#include <algorithms/TensorSlicer.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <sstream>
#include <Sisi4s.hpp>
#include <util/Parsing.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorSlicer);

template <typename F>
std::string _showVector(const std::vector<F> &v) {
  std::stringstream s;
  for (const auto &i : v) s << i << ", ";
  return "{" + s.str() + "}";
}

void TensorSlicer::run() {

  checkArgumentsOrDie({"Data", "Out", "Begin", "End"});

  std::vector<int> begin(pars::parseVector<int>(getTextArgument("Begin"))),
      end(pars::parseVector<int>(getTextArgument("End")));

  LOG(0, "TensorSlicer") << "begin: " << _showVector(begin) << std::endl;
  LOG(0, "TensorSlicer") << "end: " << _showVector(end) << std::endl;

  auto T(getTensorArgument<double>("Data"));

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

  allocatedTensorArgument<double>("Out", t);
}
