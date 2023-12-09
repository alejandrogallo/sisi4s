#include <algorithms/ApproximateParticleHoleCoulombVertex.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(ApproximateParticleHoleCoulombVertex,
        SPEC_IN({"FullParticleHoleCoulombVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"ParticleHoleCoulombVertexSingularVectors",
                 SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT({"ParticleHoleCoulombVertex",
                  SPEC_VAROUT("TODO: DOC", Tensor<complex> *)}));

IMPLEMENT_ALGORITHM(ApproximateParticleHoleCoulombVertex) {
  Tensor<complex> *GammaGai(
      in.get<Tensor<complex> *>("FullParticleHoleCoulombVertex"));
  Tensor<complex> *UGF(
      in.get<Tensor<complex> *>("ParticleHoleCoulombVertexSingularVectors"));
  int lens[] = {UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2]};
  int syms[] = {NS, NS, NS};
  Tensor<complex> *GammaFai =
      new Tensor<complex>(3, lens, syms, *GammaGai->wrld, "GammaFai");
  out.set<Tensor<complex> *>("ParticleHoleCoulombVertex", GammaFai);
  (*GammaFai)["Fai"] = (*GammaGai)["Gai"] * (*UGF)["GF"];
}

void ApproximateParticleHoleCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGai(
      in.get<DryTensor<complex> *>("FullParticleHoleCoulombVertex"));
  DryTensor<complex> *UGF(
      in.get<DryTensor<complex> *>("ParticleHoleCoulombVertexSingularVectors"));
  int lens[] = {UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2]};
  int syms[] = {NS, NS, NS};
  DryTensor<complex> *GammaFai =
      new DryTensor<complex>(3, lens, syms, SOURCE_LOCATION);
  out.set<DryTensor<complex> *>("ParticleHoleCoulombVertex", GammaFai);
}
