#include <algorithms/FockMatrixFromCoulombIntegrals.hpp>
#include <string>
#include <vector>
#include <math/MathFunctions.hpp>
#include <algorithm>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <map>
#include <set>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(FockMatrixFromCoulombIntegrals);

void FockMatrixFromCoulombIntegrals::run() {

  const bool unrestricted(getIntegerArgument("unrestricted", 0) == 1);

  LOG(0, "FockMatrixFromCoulombIntegrals")
    << "unrestricted? " << unrestricted
    << std::endl;

  const auto phph(getTensorArgument<double>("PHPHCoulombIntegrals"));
  int No(phph->lens[1]), Nv(phph->lens[0]);
  std::vector<int> vv({Nv, Nv}), oo({No, No}), ov({No, Nv}), vo({Nv, No});
  std::vector<int> syms({NS, NS});

  // ij: HH HHHH
  auto fij(new CTF::Tensor<double>(2, oo.data(), syms.data(), *Cc4s::world));
  const auto hh(getTensorArgument<double>("HHMatrix"));
  const auto hhhh(getTensorArgument<double>("HHHHCoulombIntegrals"));
  (*fij)["ij"]  = (*hh)["ij"];
  (*fij)["ij"] += (unrestricted ? 1.0 : 2.0)*(*hhhh)["ikjk"];
  (*fij)["ij"] += (-1.0)*(*hhhh)["ikkj"];
  allocatedTensorArgument<double>("HHFockMatrix", fij);

  // ab: PP PHPH PHHP
  auto fab(new CTF::Tensor<double>(2, vv.data(), syms.data(), *Cc4s::world));
  const auto pp(getTensorArgument<double>("PPMatrix"));
  const auto phhp(getTensorArgument<double>("PHHPCoulombIntegrals"));
  (*fab)["ab"] = (*pp)["ab"];
  (*fab)["ab"] += (unrestricted ? 1.0 : 2.0)*(*phph)["akbk"];
  (*fab)["ab"] += (-1.0)*(*phhp)["akkb"];
  allocatedTensorArgument<double>("PPFockMatrix", fab);

  // ai: PH PHHH
  auto fai(new CTF::Tensor<double>(2, vo.data(), syms.data(), *Cc4s::world));
  const auto ph(getTensorArgument<double>("PHMatrix"));
  const auto phhh(getTensorArgument<double>("PHHHCoulombIntegrals"));
  (*fai)["ai"] = (*ph)["ai"];
  (*fai)["ai"] += (unrestricted ? 1.0 : 2.0)*(*phhh)["akik"];
  (*fai)["ai"] += (-1.0)*(*phhh)["akki"];
  allocatedTensorArgument<double>("PHFockMatrix", fai);

  // ia: HP HHPH HHHP
  auto fia(new CTF::Tensor<double>(2, ov.data(), syms.data(), *Cc4s::world));
  (*fia)["ia"] = (*fai)["ai"];
  allocatedTensorArgument<double>("HPFockMatrix", fia);

  auto epsi(new CTF::Tensor<double>(1, oo.data(), syms.data(), *Cc4s::world));
  (*epsi)["i"] = (*fij)["ii"];
  allocatedTensorArgument<double>("HoleEigenEnergies", epsi);

  auto epsa(new CTF::Tensor<double>(1, vv.data(), syms.data(), *Cc4s::world));
  (*epsa)["a"] = (*fab)["aa"];
  allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);

  CTF::Scalar<double> energy;
  energy[""]  = (unrestricted ?  1.0 :  2.0) * (*hh)["ii"];
  energy[""] += (unrestricted ?  0.5 :  2.0) * (*hhhh)["ikik"];
  energy[""] += (unrestricted ? -0.5 : -1.0) * (*hhhh)["ikki"];
  const double dEnergy(energy.get_val());

  LOG(0, "FockMatrixFromCoulombIntegrals")
    << "energie = " << dEnergy << std::endl;

}
