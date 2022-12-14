#include <algorithms/DoublesAmplitudesFromVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(DoublesAmplitudesFromVertex);

DoublesAmplitudesFromVertex::DoublesAmplitudesFromVertex(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

DoublesAmplitudesFromVertex::~DoublesAmplitudesFromVertex() {}

void DoublesAmplitudesFromVertex::run() {
  // read the amplitudes vertex YLai
  Tensor<complex> *YLai(getTensorArgument<complex>("DoublesAmplitudesVertex"));

  // get Nv,No
  int Nv(YLai->lens[1]);
  int No(YLai->lens[2]);

  // allocate amplitudes
  int syms[] = {NS, NS, NS, NS};
  int vvoo[] = {Nv, Nv, No, No};
  Tensor<double> *Tabij(new Tensor<>(4, vvoo, syms, *Sisi4s::world, "Tabij"));

  // split YLai into real and imaginary parts
  Tensor<double> realYLai(3, YLai->lens, YLai->sym, *YLai->wrld, "realYLai");
  Tensor<double> imagYLai(3, YLai->lens, YLai->sym, *YLai->wrld, "imagYLai");
  fromComplexTensor(*YLai, realYLai, imagYLai);

  (*Tabij)["abij"] = realYLai["Lai"] * realYLai["Lbj"];
  (*Tabij)["abij"] -= imagYLai["Lai"] * imagYLai["Lbj"];

  allocatedTensorArgument("DoublesAmplitudes", Tabij);
}

void DoublesAmplitudesFromVertex::dryRun() {
  // read the Coulomb vertex GammaGqr
  DryTensor<complex> *YLai(getTensorArgument<complex, DryTensor<complex>>(
      "DoublesAmplitudesVertex"));

  // get Nv,No
  int Nv(YLai->lens[1]);
  int No(YLai->lens[2]);

  // allocate amplitudes
  int syms[] = {NS, NS, NS, NS};
  int vvoo[] = {Nv, Nv, No, No};
  DryTensor<> *Tabij(new DryTensor<>(4, vvoo, syms, SOURCE_LOCATION));
  Tabij->use();

  // split YLai into real and imaginary parts
  DryTensor<> realYLai(3,
                       YLai->lens.data(),
                       YLai->syms.data(),
                       SOURCE_LOCATION);
  DryTensor<> imagYLai(3,
                       YLai->lens.data(),
                       YLai->syms.data(),
                       SOURCE_LOCATION);
}
