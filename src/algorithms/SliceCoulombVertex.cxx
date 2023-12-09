#include <algorithms/SliceCoulombVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(SliceCoulombVertex,
        SPEC_IN({"CoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
                {"HoleEigenEnergies",
                 SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
        SPEC_OUT({"ParticleHoleCoulombVertex",
                  SPEC_VAROUT("TODO: DOC", Tensor<complex> *)}));

IMPLEMENT_ALGORITHM(SliceCoulombVertex) {
  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr(in.get<Tensor<complex> *>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Np(GammaGqr->lens[1]);

  // Allocate and compute GammaGai
  int GaiStart[] = {0, No, 0};
  int GaiEnd[] = {NG, Np, No};
  Tensor<complex> *GammaGai(
      new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd)));
  out.set<Tensor<complex> *>("ParticleHoleCoulombVertex", GammaGai);
}

void SliceCoulombVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(in.get<DryTensor<complex> *>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(in.get<DryTensor<double> *>("HoleEigenEnergies"));
  DryTensor<> *epsa(in.get<DryTensor<double> *>("ParticleEigenEnergies"));

  // Compute the No,Nv,NG
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiLens[] = {NG, Nv, No};
  int GaiSyms[] = {NS, NS, NS};
  DryTensor<complex> *GammaGai(new DryTensor<complex>(3, GaiLens, GaiSyms));

  out.set<DryTensor<complex> *>("ParticleHoleCoulombVertex", GammaGai);
}
