#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <typeinfo>
#include <vector>

#include <Sisi4s.hpp>
#include <algorithms/TensorAntisymmetrizer.hpp>
#include <math/MathFunctions.hpp>
#include <util/Integrals.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>
#include <util/Tensor.hpp>

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

static const std::type_info &check_type(Arguments &in) {
  const std::vector<IntegralInfo> infos = get_integral_infos();
  for (const auto &integral : infos)
    if (in.present(integral.name)) {
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
static void run(Arguments &in) {
  const std::vector<IntegralInfo> infos = get_integral_infos();

  // Only copy integrals that are not self-antisymmetrizable
  std::map<std::string, PTR(Tensor<F>)> integralCopies;
  for (const auto &integral : infos) {
    if (in.present(integral.name) && !isSelfAntisymmetrizable(integral)) {
      LOG(1, "TensorAntisymmetrizer")
          << "Copying " << integral.name << std::endl;
      integralCopies[integral.name] =
          NEW(Tensor<F>, in.get<Tensor<F> *>(integral.name));
    }
  }

  for (const auto &integral : infos) {
    bool antisymmetrized(false);
    if (!in.present(integral.name)) continue;

    auto antigrals(integral.getAntisymmetrizers());
    // sort antigrals so that integral be the first
    // if it is in antigrals
    sort(antigrals.begin(),
         antigrals.end(),
         [&](const IntegralInfo &a, const IntegralInfo &b) {
           return a.name == integral.name;
         });

    for (const auto &antigral : antigrals) {
      if (in.present(antigral.name)) {
        LOG(1, "TensorAntisymmetrizer")
            << integral.name << " from " << antigral.name << std::endl;
        LOG(1, "TensorAntisymmetrizer")
            << "  " << integral.name << "[" << integral.ids
            << "] -= " << antigral.name << "[" << antigral.ids << "]"
            << std::endl;
        auto inteCtf(in.get<Tensor<F> *>(integral.name));
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

using FSPEC = double;
DEFSPEC(
    TensorAntisymmetrizer,
    SPEC_IN({"mode", SPEC_ONE_OF("TODO: DOC", std::string, "up", "down")},
            {"HHHHCoulombIntegrals", SPEC_VARIN("Vijkl", Tensor<FSPEC *>())},
            {"HHHPCoulombIntegrals", SPEC_VARIN("Vijka", Tensor<FSPEC *>())},
            {"HHPHCoulombIntegrals", SPEC_VARIN("Vijak", Tensor<FSPEC *>())},
            {"HHPPCoulombIntegrals", SPEC_VARIN("Vijab", Tensor<FSPEC *>())},
            {"HPHHCoulombIntegrals", SPEC_VARIN("Viajk", Tensor<FSPEC *>())},
            {"HPHPCoulombIntegrals", SPEC_VARIN("Viajb", Tensor<FSPEC *>())},
            {"HPPHCoulombIntegrals", SPEC_VARIN("Viabj", Tensor<FSPEC *>())},
            {"HPPPCoulombIntegrals", SPEC_VARIN("Viabc", Tensor<FSPEC *>())},
            {"PHHHCoulombIntegrals", SPEC_VARIN("Vaijk", Tensor<FSPEC *>())},
            {"PHHPCoulombIntegrals", SPEC_VARIN("Vaijb", Tensor<FSPEC *>())},
            {"PHPHCoulombIntegrals", SPEC_VARIN("Vaibj", Tensor<FSPEC *>())},
            {"PHPPCoulombIntegrals", SPEC_VARIN("Vaibc", Tensor<FSPEC *>())},
            {"PPHHCoulombIntegrals", SPEC_VARIN("Vabij", Tensor<FSPEC *>())},
            {"PPHPCoulombIntegrals", SPEC_VARIN("Vabic", Tensor<FSPEC *>())},
            {"PPPHCoulombIntegrals", SPEC_VARIN("Vabci", Tensor<FSPEC *>())},
            {"PPPPCoulombIntegrals", SPEC_VARIN("Vabcd", Tensor<FSPEC *>())},
            {"left", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},
            {"right", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)}),
    SPEC_OUT());

template <typename F>
static void run2(Arguments &in) {
  auto &left = *in.get<Tensor<F> *>("left"),
       &right = *in.get<Tensor<F> *>("right");

  const std::string mode = in.get<std::string>("mode");

  if (left.order != 4)
    throw "TensorAntisymmetrizer2 can only antisymmetrize 4 index tensors";

  if (mode == "up") left["abij"] -= right["baij"];
  else left["abij"] -= right["abji"];
}

IMPLEMENT_ALGORITHM(TensorAntisymmetrizer) {
  if (in.present("left") && in.present("right")) {
    if (in.is_of_type<Tensor<double> *>("right")
        && in.is_of_type<Tensor<double> *>("left")) {
      ::run2<double>(in);
    } else {
      ::run2<sisi4s::complex>(in);
    }
  } else {
    if (typeid(double) == check_type(in)) {
      LOG(1, "TensorAntisymmetrizer")
          << "real integrals detected " << std::endl;
      ::run<double>(in);
    } else {
      LOG(1, "TensorAntisymmetrizer")
          << "complex integrals detected " << std::endl;
      ::run<sisi4s::complex>(in);
    }
  }
}
