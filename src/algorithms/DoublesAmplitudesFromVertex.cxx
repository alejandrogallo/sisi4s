#include <algorithms/DoublesAmplitudesFromVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(DoublesAmplitudesFromVertex,
        SPEC_IN({"DoublesAmplitudesVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
        SPEC_OUT({"DoublesAmplitudes",
                  SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(DoublesAmplitudesFromVertex) {
  // read the amplitudes vertex YLai
  Tensor<sisi4s::complex> *YLai(
      in.get<Tensor<sisi4s::complex> *>("DoublesAmplitudesVertex"));

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

  out.set<Tensor<double> *>("DoublesAmplitudes", Tabij);
}

void DoublesAmplitudesFromVertex::dryRun() {
  // read the Coulomb vertex GammaGqr
  DryTensor<sisi4s::complex> *YLai(
      in.get<DryTensor<sisi4s::complex> *>("DoublesAmplitudesVertex"));

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
