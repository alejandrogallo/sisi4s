#include <algorithms/SingleParticleOccupancies.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;


DEFSPEC(
    SingleParticleOccupancies,
    SPEC_IN({"DoublesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT({"HoleOccupancies", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"ParticleOccupancies",
              SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(SingleParticleOccupancies) {
  // Read the DRCCD amplitudes Tabij
  Tensor<double> *Tabij(in.get<Tensor<double> *>("DoublesAmplitudes"));

  int no(Tabij->lens[2]), nv(Tabij->lens[0]);
  // create particle and hole occupancies
  Vector<> *Ni(new Vector<>(no, *Tabij->wrld, "Ni"));
  Vector<> *Na(new Vector<>(nv, *Tabij->wrld, "No"));

  // calculate <T|T> which is <Psi|Psi> in a linearized CC theory
  Scalar<> TT(*Tabij->wrld);
  TT.set_name("TT");
  TT[""] = (*Tabij)["abij"] * (*Tabij)["abij"];

  // calculate <Psi|Na|Psi>
  (*Na)["a"] = +2.0 * (*Tabij)["abij"] * (*Tabij)["abij"];
  // calculate <Psi|Na|Psi> / <Psi|Psi>
  Bivar_Function<> fDivide(&divide<double>);
  Tensor<double> Da(false, *Na);
  Da["a"] = TT[""];
  Na->contract(1.0, *Na, "a", Da, "a", 0.0, "a", fDivide);

  // calculate <Psi|Ni|Psi>
  (*Ni)["i"] = -2.0 * (*Tabij)["abij"] * (*Tabij)["abij"];
  // calculate <Psi|Na|Psi> / <Psi|Psi>
  Tensor<double> Di(false, *Ni);
  Di["i"] = TT[""];
  Ni->contract(1.0, *Ni, "i", Di, "i", 0.0, "i", fDivide);
  (*Ni)["i"] += 1.0;

  out.set<Tensor<double> *>("ParticleOccupancies", Na);
  out.set<Tensor<double> *>("HoleOccupancies", Ni);
}

void SingleParticleOccupancies::dryRun() {
  // Read the DRCCD amplitudes Tabij
  DryTensor<> *Tabij(in.get<DryTensor<double> *>("DoublesAmplitudes"));

  int no(Tabij->lens[2]), nv(Tabij->lens[0]);
  // create particle and hole occupancies
  DryVector<> *Ni(new DryVector<>(no));
  DryVector<> *Na(new DryVector<>(nv));

  // calculate <T|T> which is <Psi|Psi> in a linearized CC theory
  DryScalar<> TT();

  // calculate <Psi|Na|Psi>
  DryTensor<> Da(*Na);
  DryTensor<> Di(*Ni);

  out.set<DryTensor<double> *>("ParticleOccupancies", Na);
  out.set<DryTensor<double> *>("HoleOccupancies", Ni);
}
