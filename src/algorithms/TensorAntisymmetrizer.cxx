#include <Sisi4s.hpp>
#include <algorithms/TensorAntisymmetrizer.hpp>
#include <math/MathFunctions.hpp>
#include <util/Integrals.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>
#include <util/Tensor.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>

using namespace sisi4s;

static std::vector<IntegralInfo> get_integral_infos() {
  return {{"HHHHCoulombIntegrals", {NO, NO, NO, NO}, "ijkl"},
          {"HHHPCoulombIntegrals", {NO, NO, NO, NV}, "ijka"},
          {"HHPHCoulombIntegrals", {NO, NO, NV, NO}, "ijak"},
          {"HHPPCoulombIntegrals", {NO, NO, NV, NV}, "ijab"},
          {"HPHHCoulombIntegrals", {NO, NV, NO, NO}, "iajk"},
          {"HPHPCoulombIntegrals", {NO, NV, NO, NV}, "iajb"},
          {"HPPHCoulombIntegrals", {NO, NV, NV, NO}, "iabj"},
          {"HPPPCoulombIntegrals", {NO, NV, NV, NV}, "iabc"},
          {"PHHHCoulombIntegrals", {NV, NO, NO, NO}, "aijk"},
          {"PHHPCoulombIntegrals", {NV, NO, NO, NV}, "aijb"},
          {"PHPHCoulombIntegrals", {NV, NO, NV, NO}, "aibj"},
          {"PHPPCoulombIntegrals", {NV, NO, NV, NV}, "aibc"},
          {"PPHHCoulombIntegrals", {NV, NV, NO, NO}, "abij"},
          {"PPHPCoulombIntegrals", {NV, NV, NO, NV}, "abic"},
          {"PPPHCoulombIntegrals", {NV, NV, NV, NO}, "abci"},
          {"PPPPCoulombIntegrals", {NV, NV, NV, NV}, "abcd"}};
}

static const std::type_info &check_type(Algorithm &alg) {
  const std::vector<IntegralInfo> infos = get_integral_infos();
  for (const auto &integral : infos)
    if (alg.in.present(integral.name)) {
      bool real_tensor_data = in.is_of_type<Tensor<double> *>(integral.name);
      if (real_tensor_data) return typeid(double);
      bool imag_tensor_data =
          in.is_of_type<Tensor<sisi4s::complex>>(integral.name);
      if (imag_tensor_data) return typeid(sisi4s::complex);
    }
  throw "Could not detect type of integrals in TensorAntisymmetrizer";
}

static bool isSelfAntisymmetrizable(const IntegralInfo &i) {
  for (const auto &j : i.getAntisymmetrizers()) {
    if (j.name == i.name) return true;
  }
  return false;
}

template <typename F>
static void run(Algorithm &alg) {
  const std::vector<IntegralInfo> infos = get_integral_infos();

  // Only copy integrals that are not self-antisymmetrizable
  std::map<std::string, PTR(Tensor<F>)> integralCopies;
  for (const auto &integral : infos) {
    if (alg.isArgumentGiven(integral.name)
        && !isSelfAntisymmetrizable(integral)) {
      LOG(1, "TensorAntisymmetrizer")
          << "Copying " << integral.name << std::endl;
      integralCopies[integral.name] =
          NEW(Tensor<F>, *alg.in.get<Tensor<F> *>(integral.name));
    }
  }

  for (const auto &integral : infos) {
    bool antisymmetrized(false);
    if (!alg.isArgumentGiven(integral.name)) continue;

    auto antigrals(integral.getAntisymmetrizers());
    // sort antigrals so that integral be the first
    // if it is in antigrals
    sort(antigrals.begin(),
         antigrals.end(),
         [&](const IntegralInfo &a, const IntegralInfo &b) {
           return a.name == integral.name;
         });

    for (const auto &antigral : antigrals) {
      if (alg.isArgumentGiven(antigral.name)) {
        LOG(1, "TensorAntisymmetrizer")
            << integral.name << " from " << antigral.name << std::endl;
        LOG(1, "TensorAntisymmetrizer")
            << "  " << integral.name << "[" << integral.ids
            << "] -= " << antigral.name << "[" << antigral.ids << "]"
            << std::endl;
        auto inteCtf(alg.in.get<Tensor<F> *>(integral.name));
        Tensor<F> *antiCtf;
        if (antigral.name == integral.name) {
          antiCtf = inteCtf;
        } else {
          antiCtf = integralCopies[antigral.name].get();
        }
        (*inteCtf)[integral.ids.data()] -= (*antiCtf)[antigral.ids.data()];
        antisymmetrized = true;
        // go to the next integral
        break;
      }
    }

    // if integral.name is not antisymmetrized, then something went
    // wrong, we don't have enough integrals to do this, panic!
    if (!antisymmetrized) {
      LOG(1, "TensorAntisymmetrizer")
          << "error: " << integral.name << " could not be antisymmetrized"
          << std::endl;
      LOG(1, "TensorAntisymmetrizer") << "error: possibilities: " << std::endl;
      for (const auto &antigral : integral.getAntisymmetrizers()) {
        LOG(1, "TensorAntisymmetrizer") << "> " << antigral.name << std::endl;
      }
      throw new EXCEPTION("Antisymmetrization error");
    }
  }
}


DEFSPEC(TensorAntisymmetrizer,
        SPEC_IN({"mode", SPEC_VALUE_DEF("TODO: DOC", std::string, "up")},
                {integral.name, SPEC_VARIN("TODO: DOC", Tensor<F> *)},
                {"left", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
                {"right", SPEC_VARIN("TODO: DOC", Tensor<F> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(TensorAntisymmetrizer) {
  if (typeid(double) == check_type(*this)) {
    LOG(1, "TensorAntisymmetrizer") << "real integrals detected " << std::endl;
    ::run<double>(*this);
  } else {
    LOG(1, "TensorAntisymmetrizer")
        << "complex integrals detected " << std::endl;
    ::run<sisi4s::complex>(*this);
  }
}

template <typename F>
static void run2(Algorithm &alg) {
  auto &left = *alg.in.get<Tensor<F> *>("left"),
       &right = *alg.in.get<Tensor<F> *>("right");

  const std::string mode = alg.in.get<std::string>("mode", "up");

  if (left.order != 4)
    throw "TensorAntisymmetrizer2 can only antisymmetrize 4 index tensors";

  if (mode == "up") left["abij"] -= right["baij"];
  else left["abij"] -= right["abji"];
}


DEFSPEC(TensorAntisymmetrizer,
        SPEC_IN({"mode", SPEC_VALUE_DEF("TODO: DOC", std::string, "up")},
                {integral.name, SPEC_VARIN("TODO: DOC", Tensor<F> *)},
                {"left", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
                {"right", SPEC_VARIN("TODO: DOC", Tensor<F> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(TensorAntisymmetrizer2) {
  if (in.is_of_type<Tensor<double> *>("Data")) {
    LOG(1, "TensorAntisymmetrizer2") << "real integrals detected " << std::endl;
    ::run2<double>(*this);
  } else {
    LOG(1, "TensorAntisymmetrizer2")
        << "complex integrals detected " << std::endl;
    ::run2<sisi4s::complex>(*this);
  }
}
