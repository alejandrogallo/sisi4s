#include <algorithms/CcsdtEnergyFromCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <array>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdtEnergyFromCoulombIntegrals);


template <typename F>
PTR(FockVector<F>)
CcsdtEnergyFromCoulombIntegrals::getResiduumTemplate(const int iterationStep,
                                                     const PTR(const FockVector<F>)
                                                     &amplitudes) {

  auto
    *epsi = getTensorArgument<double, Tensor<double> >("HoleEigenEnergies"),
    *epsa = getTensorArgument<double, Tensor<double> >("ParticleEigenEnergies");

  auto
    Vhhhh(getTensorArgument<F, Tensor<F> >("HHHHCoulombIntegrals")),
    Vpppp(getTensorArgument<F, Tensor<F> >("PPPPCoulombIntegrals")),
    Vhhhp(getTensorArgument<F, Tensor<F> >("HHHPCoulombIntegrals")),
    Vhhpp(getTensorArgument<F, Tensor<F> >("HHPPCoulombIntegrals")),
    Vhphh(getTensorArgument<F, Tensor<F> >("HPHHCoulombIntegrals")),
    Vhphp(getTensorArgument<F, Tensor<F> >("HPHPCoulombIntegrals")),
    Vhppp(getTensorArgument<F, Tensor<F> >("HPPPCoulombIntegrals")),
    Vpphh(getTensorArgument<F, Tensor<F> >("PPHHCoulombIntegrals")),
    Vpphp(getTensorArgument<F, Tensor<F> >("PPHPCoulombIntegrals")),
    Vhpph(getTensorArgument<F, Tensor<F> >("HPPHCoulombIntegrals")),
    Vphpp(getTensorArgument<F, Tensor<F> >("PHPPCoulombIntegrals")),
    Vhhph(getTensorArgument<F, Tensor<F> >("HHPHCoulombIntegrals")),
    Vppph(getTensorArgument<F, Tensor<F> >("PPPHCoulombIntegrals")),
    Vphph(getTensorArgument<F, Tensor<F> >("PHPHCoulombIntegrals")),
    Vphhp(getTensorArgument<F, Tensor<F> >("PHHPCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};

  Tensor<F>
    *Fab = new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab"),
    *Fij = new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij");


  { // build Fab and Fij
    const auto
      to_F = CTF::Transform<double, F>([](double eps, F &f) {f = eps;});
    to_F((*epsi)["i"],
         (*Fij)["ii"]);
    to_F((*epsa)["a"],
         (*Fab)["aa"]);
  }


  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  auto Tph(amplitudes->get(0));
  Tph->set_name("Tph");
  auto Tpphh(amplitudes->get(1));
  Tpphh->set_name("Tpphh");
  auto Tppphhh(amplitudes->get(2));
  Tppphhh->set_name("Tppphhh");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rph(residuum->get(0));
  Rph->set_name("Rph");
  auto Rpphh(residuum->get(1));
  Rpphh->set_name("Rpphh");
  auto Rppphhh(residuum->get(2));
  Rppphhh->set_name("Rppphhh");

  if ((iterationStep == 0) && !isArgumentGiven("initialDoublesAmplitudes")){
    LOG(1, getAbbreviation())
      << "Set initial Rpphh amplitudes to Vijab" << std::endl;
    (*Rpphh)["abij"] = (*Vhhpp)["ijab"];
    return residuum;
  }


  // TODO: insert equations


  // TODO: delete Fij and so on
  return residuum;
}

PTR(FockVector<double>)
CcsdtEnergyFromCoulombIntegrals::getResiduum(const int iterationStep,
                                             const PTR(const FockVector<double>)
                                             &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

PTR(FockVector<sisi4s::complex>)
CcsdtEnergyFromCoulombIntegrals::getResiduum(const int iterationStep,
                                             const PTR(const FockVector<sisi4s::complex>)
                                             &amplitudes) {
  return getResiduumTemplate<sisi4s::complex>(iterationStep, amplitudes);
}

