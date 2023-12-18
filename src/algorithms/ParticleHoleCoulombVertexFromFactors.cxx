#include <algorithms/ParticleHoleCoulombVertexFromFactors.hpp>
#include <math/Complex.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

#include <vector>
#include <string>
#include <memory>

using namespace sisi4s;
using namespace tcc;
using std::make_shared;

DEFSPEC(ParticleHoleCoulombVertexFromFactors,
        SPEC_IN({"CoulombFactors", SPEC_VALUE("TODO: DOC", T *)},
                {"HoleFactorOrbitals", SPEC_VALUE("TODO: DOC", T *)},
                {"ParticleFactorOrbitals", SPEC_VALUE("TODO: DOC", T *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ParticleHoleCoulombVertexFromFactors) {
  run<Tensor<sisi4s::complex>, CtfMachineTensor<sisi4s::complex>>(false);
}

void ParticleHoleCoulombVertexFromFactors::dryRun() {
  run<DryTensor<sisi4s::complex>, DryMachineTensor<sisi4s::complex>>(true);
}

// TMT is either CtfMachineTensor or DryMachineTensor
template <typename T, typename MT>
void ParticleHoleCoulombVertexFromFactors::run(const bool dryRun) {
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<complex>::create(machineTensorFactory));

  // Read the Coulomb vertex GammaGqr
  T *ctfPiiR(in.get<T *>("HoleFactorOrbitals"));
  T *ctfPiaR(in.get<T *>("ParticleFactorOrbitals"));
  T *ctfLambdaFR(in.get<T *>("CoulombFactors"));

  // for now: create tcc::Tensors from them
  // later there will only be tcc::Tensors objects stored in sisi4s
  auto PiiR(tcc->createTensor(MT::create(*ctfPiiR)));
  auto PiaR(tcc->createTensor(MT::create(*ctfPiaR)));
  auto LambdaFR(tcc->createTensor(MT::create(*ctfLambdaFR)));

  // allocate tcc::Tensor for final result
  int NF(LambdaFR->lens[0]);
  int Nv(PiaR->lens[0]);
  int No(PiiR->lens[0]);
  auto GammaFai(tcc->createTensor(std::vector<int>({NF, Nv, No}), "GammaFai"));

  // compile
  auto operation(
      tcc->compile((*GammaFai)["Fai"] <<=
                   (*LambdaFR)["FR"] * (*PiaR)["aR"] * (*PiiR)["iR"]));
  // and execute
  operation->execute();

  // for now: duplicate result
  // later Gamma will already be the object stored in sisi4s
  out.set<T *>("ParticleHoleCoulombVertex",
               new T(GammaFai->template getMachineTensor<MT>()->tensor));
}
