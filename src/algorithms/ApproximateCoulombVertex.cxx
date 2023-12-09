#include <algorithms/ApproximateCoulombVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(
    ApproximateCoulombVertex,
    SPEC_IN({"CoulombVertexSingularVectors",
             SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
            {"FullCoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
    SPEC_OUT({"CoulombVertex", SPEC_VAROUT("TODO: DOC", Tensor<complex> *)}));

IMPLEMENT_ALGORITHM(ApproximateCoulombVertex) {
  Tensor<complex> *GammaGqr(in.get<Tensor<complex> *>("FullCoulombVertex"));
  GammaGqr->set_name("GammaGqr");
  Tensor<complex> *UGF(
      in.get<Tensor<complex> *>("CoulombVertexSingularVectors"));
  Tensor<complex> UTGF(*UGF);
  conjugate(UTGF);
  int lens[] = {static_cast<int>(UGF->lens[1]),
                static_cast<int>(GammaGqr->lens[1]),
                static_cast<int>(GammaGqr->lens[2])};
  int syms[] = {NS, NS, NS};
  Tensor<complex> *GammaFqr =
      new Tensor<complex>(3, lens, syms, *GammaGqr->wrld, "GammaFqr");
  out.set<Tensor<complex> *>("CoulombVertex", GammaFqr);
  LOG(1, "ApproximateCoulombVertex")
      << "Using NF=" << UGF->lens[1] << " field variables to build "
      << GammaFqr->get_name() << " from " << GammaGqr->get_name()
      << " with NG=" << UGF->lens[0] << " grid points" << std::endl;
  (*GammaFqr)["Fqr"] = (*GammaGqr)["Gqr"] * UTGF["GF"];
}

void ApproximateCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGqr(
      in.get<DryTensor<complex> *>("FullCoulombVertex"));
  DryTensor<complex> *UGF(
      in.get<DryTensor<complex> *>("CoulombVertexSingularVectors"));
  DryTensor<complex> UTGF(*UGF);
  int lens[] = {static_cast<int>(UGF->lens[1]),
                static_cast<int>(GammaGqr->lens[1]),
                static_cast<int>(GammaGqr->lens[2])};
  int syms[] = {NS, NS, NS};
  DryTensor<complex> *GammaFqr =
      new DryTensor<complex>(3, lens, syms, SOURCE_LOCATION);
  out.set<DryTensor<complex> *>("CoulombVertex", GammaFqr);
}
