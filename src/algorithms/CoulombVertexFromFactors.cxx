#include <algorithms/CoulombVertexFromFactors.hpp>
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

DEFSPEC(CoulombVertexFromFactors,
        SPEC_IN({"CoulombFactors", SPEC_VALUE("TODO: DOC", T *)},
                {"FactorOrbitals", SPEC_VALUE("TODO: DOC", T *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(CoulombVertexFromFactors) {
  run<Tensor<sisi4s::complex>, CtfMachineTensor<sisi4s::complex>>(false);
}

void CoulombVertexFromFactors::dryRun() {
  run<DryTensor<sisi4s::complex>, DryMachineTensor<sisi4s::complex>>(true);
}

// TMT is either CtfMachineTensor or DryMachineTensor
template <typename T, typename MT>
void CoulombVertexFromFactors::run(const bool dryRun) {
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<complex>::create(machineTensorFactory));

  // Read the Coulomb vertex GammaGqr
  T *ctfPirR(in.get<T *>("FactorOrbitals"));
  T *ctfLambdaFR(in.get<T *>("CoulombFactors"));

  // for now: create tcc::Tensors from them
  // later there will only be tcc::Tensors objects stored in sisi4s
  auto PirR(tcc->createTensor(MT::create(*ctfPirR)));
  auto LambdaFR(tcc->createTensor(MT::create(*ctfLambdaFR)));

  // allocate tcc::Tensor for final result
  int NF(LambdaFR->lens[0]);
  int Np(PirR->lens[0]);
  auto GammaFqr(tcc->createTensor(std::vector<int>({NF, Np, Np}), "Gamma"));

  // compile
  auto operation(
      tcc->compile((*GammaFqr)["Fqr"] <<=
                   (*LambdaFR)["FR"] * (*PirR)["qR"] * (*PirR)["rR"]));
  // and execute
  operation->execute();

  // for now: duplicate result
  // later Gamma will already be the object stored in sisi4s
  out.set<T *>("CoulombVertex",
               new T(GammaFqr->template getMachineTensor<MT>()->tensor));
}
