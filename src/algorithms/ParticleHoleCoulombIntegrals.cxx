#include <algorithms/ParticleHoleCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(ParticleHoleCoulombIntegrals,
        SPEC_IN({"antisymmetrize", SPEC_VALUE_DEF("TODO: DOC", int64_t, 0)},
                {"ParticleHoleCoulombVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT({"HHPPCoulombIntegrals",
                  SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
                 {"PPHHCoulombIntegrals",
                  SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(ParticleHoleCoulombIntegrals) {
  // read coulomb vertex GammaGai
  Tensor<complex> *GammaGai(
      in.get<Tensor<complex> *>("ParticleHoleCoulombVertex"));

  // allocate real and imag part of GammaGai
  Tensor<double> realGammaGai(3,
                              GammaGai->lens,
                              GammaGai->sym,
                              *GammaGai->wrld,
                              "RealGammaGai");
  Tensor<double> imagGammaGai(3,
                              GammaGai->lens,
                              GammaGai->sym,
                              *GammaGai->wrld,
                              "ImagGammaGai");

  // split into real and imaginary parts
  fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

  // allocate coulomb integrals
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
  int vvoo[] = {Nv, Nv, No, No};
  int oovv[] = {No, No, Nv, Nv};
  int syms[] = {NS, NS, NS, NS};
  Tensor<double> *Vabij(
      isArgumentGiven("PPHHCoulombIntegrals")
          ? new Tensor<double>(4, vvoo, syms, *Sisi4s::world, "Vabij")
          : nullptr);
  Tensor<double> *Vijab(
      isArgumentGiven("HHPPCoulombIntegrals")
          ? new Tensor<double>(4, oovv, syms, *Sisi4s::world, "Vijab")
          : nullptr);
  if (Vabij) {
    out.set<Tensor<double> *>("PPHHCoulombIntegrals", Vabij);
    (*Vabij)["abij"] = realGammaGai["gai"] * realGammaGai["gbj"];
    (*Vabij)["abij"] += imagGammaGai["gai"] * imagGammaGai["gbj"];
  }
  if (Vijab) {
    // FIXME: allow Complex integrals
    out.set<Tensor<double> *>("HHPPCoulombIntegrals", Vijab);
    (*Vijab)["ijab"] = (*Vabij)["abij"];
  }

  int antisymmetrize(in.get<int64_t>("antisymmetrize", 0));
  if (antisymmetrize) {
    LOG(0, "CoulombIntegrals")
        << "Calculating antisymmetrized integrals" << std::endl;
    if (Vabij) (*Vabij)["abij"] -= (*Vabij)["abji"];
    if (Vijab) (*Vijab)["ijab"] -= (*Vijab)["jiab"];
  }
}

void ParticleHoleCoulombIntegrals::dryRun() {
  DryTensor<complex> *GammaGai(
      in.get<DryTensor<complex> *>("ParticleHoleCoulombVertex"));

  // Compute the No,Nv,NG,Np
  int NG(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = {NS, NS, NS, NS};
  int vvoo[] = {Nv, Nv, No, No};

  DryTensor<> *Vabij(new DryTensor<>(4, vvoo, syms));

  out.set<Tensor<double> *>("PPHHCoulombIntegrals", Vabij);

  // Allocate and realGammaGai and imagGammaGai
  int GaiLens[] = {NG, Nv, No};

  DryTensor<> realGammaGai(3, GaiLens, syms);
  DryTensor<> imagGammaGai(3, GaiLens, syms);
}
