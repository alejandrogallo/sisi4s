#include <algorithms/SliceCoulombVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(SliceCoulombVertex);

void SliceCoulombVertex::run() {
  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<double> *epsi(getTensorArgument<>("HoleEigenEnergies"));

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Np(GammaGqr->lens[1]);

  // Allocate and compute GammaGai
  int GaiStart[] = {0, No, 0};
  int GaiEnd[] = {NG, Np, No};
  Tensor<complex> *GammaGai(
      new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd)));
  allocatedTensorArgument<complex>("ParticleHoleCoulombVertex", GammaGai);
}

void SliceCoulombVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(
      getTensorArgument<complex, DryTensor<complex>>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
      getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(
      getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute the No,Nv,NG
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiLens[] = {NG, Nv, No};
  int GaiSyms[] = {NS, NS, NS};
  DryTensor<complex> *GammaGai(new DryTensor<complex>(3, GaiLens, GaiSyms));

  allocatedTensorArgument<complex, DryTensor<complex>>(
      "ParticleHoleCoulombVertex",
      GammaGai);
}
