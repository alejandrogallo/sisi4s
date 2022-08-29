#include <algorithms/ParticleHoleCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <util/CTF.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombIntegrals);

ParticleHoleCoulombIntegrals::ParticleHoleCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ParticleHoleCoulombIntegrals::~ParticleHoleCoulombIntegrals() {
}

void ParticleHoleCoulombIntegrals::run() {
  // read coulomb vertex GammaGai
  Tensor<complex> *GammaGai(getTensorArgument<complex>
                            ("ParticleHoleCoulombVertex"));

  // allocate real and imag part of GammaGai
  Tensor<> realGammaGai(3, GammaGai->lens, GammaGai->sym, 
                        *GammaGai->wrld, "RealGammaGai");
  Tensor<> imagGammaGai(3, GammaGai->lens, GammaGai->sym, 
                        *GammaGai->wrld, "ImagGammaGai");

  // split into real and imaginary parts
  fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

  // allocate coulomb integrals
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
  int vvoo[] = { Nv, Nv, No, No };
  int oovv[] = { No, No, Nv, Nv };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Vabij(
    isArgumentGiven("PPHHCoulombIntegrals") ?
    new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij") : nullptr
  );
  Tensor<> *Vijab(
    isArgumentGiven("HHPPCoulombIntegrals") ?
    new Tensor<>(4, oovv, syms, *Cc4s::world, "Vijab") : nullptr
  );
  if (Vabij) {
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
    (*Vabij)["abij"] =  realGammaGai["gai"] * realGammaGai["gbj"];
    (*Vabij)["abij"] += imagGammaGai["gai"] * imagGammaGai["gbj"];
  }
  if (Vijab) {
    // FIXME: allow Complex integrals
    allocatedTensorArgument("HHPPCoulombIntegrals", Vijab);
    (*Vijab)["ijab"] = (*Vabij)["abij"];
  }

  int antisymmetrize(getIntegerArgument("antisymmetrize", 0));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals") << "Calculating antisymmetrized integrals"
      << std::endl;
    if (Vabij) (*Vabij)["abij"] -= (*Vabij)["abji"];
    if (Vijab) (*Vijab)["ijab"] -= (*Vijab)["jiab"];
  }
}

void ParticleHoleCoulombIntegrals::dryRun() {
  DryTensor<complex> *GammaGai(getTensorArgument<complex, 
                               DryTensor<complex>>
                               ("ParticleHoleCoulombVertex"));

  // Compute the No,Nv,NG,Np
  int NG(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };

  DryTensor<> *Vabij(new DryTensor<>(4, vvoo, syms));

  allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);

  // Allocate and realGammaGai and imagGammaGai
  int GaiLens[]   = {NG,Nv,No};

  DryTensor<> realGammaGai(3, GaiLens, syms);
  DryTensor<> imagGammaGai(3, GaiLens, syms);
}
