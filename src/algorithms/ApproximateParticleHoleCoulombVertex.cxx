#include <algorithms/ApproximateParticleHoleCoulombVertex.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(ApproximateParticleHoleCoulombVertex,
        SPEC_IN({"FullParticleHoleCoulombVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},
                {"ParticleHoleCoulombVertexSingularVectors",
                 SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
        SPEC_OUT({"ParticleHoleCoulombVertex",
                  SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)}));

IMPLEMENT_ALGORITHM(ApproximateParticleHoleCoulombVertex) {
  Tensor<sisi4s::complex> *GammaGai(
      in.get<Tensor<sisi4s::complex> *>("FullParticleHoleCoulombVertex"));
  Tensor<sisi4s::complex> *UGF(in.get<Tensor<sisi4s::complex> *>(
      "ParticleHoleCoulombVertexSingularVectors"));
  int lens[] = {UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2]};
  int syms[] = {NS, NS, NS};
  Tensor<sisi4s::complex> *GammaFai =
      new Tensor<sisi4s::complex>(3, lens, syms, *GammaGai->wrld, "GammaFai");
  out.set<Tensor<sisi4s::complex> *>("ParticleHoleCoulombVertex", GammaFai);
  (*GammaFai)["Fai"] = (*GammaGai)["Gai"] * (*UGF)["GF"];
}

void ApproximateParticleHoleCoulombVertex::dryRun() {
  DryTensor<sisi4s::complex> *GammaGai(
      in.get<DryTensor<sisi4s::complex> *>("FullParticleHoleCoulombVertex"));
  DryTensor<sisi4s::complex> *UGF(in.get<DryTensor<sisi4s::complex> *>(
      "ParticleHoleCoulombVertexSingularVectors"));
  int lens[] = {UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2]};
  int syms[] = {NS, NS, NS};
  DryTensor<sisi4s::complex> *GammaFai =
      new DryTensor<sisi4s::complex>(3, lens, syms, SOURCE_LOCATION);
  out.set<DryTensor<sisi4s::complex> *>("ParticleHoleCoulombVertex", GammaFai);
}
