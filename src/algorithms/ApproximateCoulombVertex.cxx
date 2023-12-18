#include <algorithms/ApproximateCoulombVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(ApproximateCoulombVertex,
        SPEC_IN({"CoulombVertexSingularVectors",
                 SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},
                {"FullCoulombVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
        SPEC_OUT({"CoulombVertex",
                  SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)}));

IMPLEMENT_ALGORITHM(ApproximateCoulombVertex) {
  Tensor<sisi4s::complex> *GammaGqr(
      in.get<Tensor<sisi4s::complex> *>("FullCoulombVertex"));
  GammaGqr->set_name("GammaGqr");
  Tensor<sisi4s::complex> *UGF(
      in.get<Tensor<sisi4s::complex> *>("CoulombVertexSingularVectors"));
  Tensor<sisi4s::complex> UTGF(*UGF);
  conjugate(UTGF);
  int lens[] = {static_cast<int>(UGF->lens[1]),
                static_cast<int>(GammaGqr->lens[1]),
                static_cast<int>(GammaGqr->lens[2])};
  int syms[] = {NS, NS, NS};
  Tensor<sisi4s::complex> *GammaFqr =
      new Tensor<sisi4s::complex>(3, lens, syms, *GammaGqr->wrld, "GammaFqr");
  out.set<Tensor<sisi4s::complex> *>("CoulombVertex", GammaFqr);
  LOG(1, "ApproximateCoulombVertex")
      << "Using NF=" << UGF->lens[1] << " field variables to build "
      << GammaFqr->get_name() << " from " << GammaGqr->get_name()
      << " with NG=" << UGF->lens[0] << " grid points" << std::endl;
  (*GammaFqr)["Fqr"] = (*GammaGqr)["Gqr"] * UTGF["GF"];
}

void ApproximateCoulombVertex::dryRun() {
  DryTensor<sisi4s::complex> *GammaGqr(
      in.get<DryTensor<sisi4s::complex> *>("FullCoulombVertex"));
  DryTensor<sisi4s::complex> *UGF(
      in.get<DryTensor<sisi4s::complex> *>("CoulombVertexSingularVectors"));
  DryTensor<sisi4s::complex> UTGF(*UGF);
  int lens[] = {static_cast<int>(UGF->lens[1]),
                static_cast<int>(GammaGqr->lens[1]),
                static_cast<int>(GammaGqr->lens[2])};
  int syms[] = {NS, NS, NS};
  DryTensor<sisi4s::complex> *GammaFqr =
      new DryTensor<sisi4s::complex>(3, lens, syms, SOURCE_LOCATION);
  out.set<DryTensor<sisi4s::complex> *>("CoulombVertex", GammaFqr);
}
