#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/MultiCombinations.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalMp2EnergyFromCoulombIntegrals);

ThermalMp2EnergyFromCoulombIntegrals::ThermalMp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalMp2EnergyFromCoulombIntegrals::~ThermalMp2EnergyFromCoulombIntegrals() {
}

void ThermalMp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  beta = 1/getRealArgument("Temperature");

  // compute \Delta^{ab}_{ij} = eps_a + eps_b - eps_i - eps_j
  Dabij = NEW(Tensor<>, false, *Vabij);
  (*Dabij)["abij"] =  (*epsa)["a"];
  (*Dabij)["abij"] += (*epsa)["b"];
  (*Dabij)["abij"] -= (*epsi)["i"];
  (*Dabij)["abij"] -= (*epsi)["j"];

  computeFreeEnergy();
  computeEnergyMoments();
}

void ThermalMp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("ThermalPPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}

void ThermalMp2EnergyFromCoulombIntegrals::computeFreeEnergy() {
  // MP2 contribution (2nd order FT-MBPT)
  Tensor<> Tabij(false, *getTensorArgument("ThermalPPHHCoulombIntegrals"));
  computeThermalMp2Amplitudes(Tabij, 0);
  real F(evaluateMp2("F^MP2", Tabij, -getRealArgument("Temperature")));
  // Hartree--Fock energy correction (1st order FT-MBPT)
  Tensor<> Tij(2, &Tabij.lens[2]);
  computeThermalHfAmplitudes(Tij, 0);
  allocatedTensorArgument<>("ThermalHfAmplitudes", new Tensor<>(Tij));
  F += evaluateHf("F^HF", Tij, -getRealArgument("Temperature"));
  // TODO: H_0 contribution: log(Z_0(beta)) = sum_1^Np log(1+exp(-eps_p*beta))
  setRealArgument("ThermalMp2FreeEnergy", F);
}

void ThermalMp2EnergyFromCoulombIntegrals::computeEnergyMoments() {
  unsigned int n( std::max(getIntegerArgument("maxEnergyMoment", 2),1l) );
  std::vector<real> energyMoments(n);

  for (unsigned int k(1); k <= n; ++k) {
    // MP2 contribution (2nd FT-MBPT)
    std::stringstream mp2Contribution;
    mp2Contribution << k << ".central moment of H^MP2";
    Tensor<> Tabij(false, *getTensorArgument("ThermalPPHHCoulombIntegrals"));
    computeThermalMp2Amplitudes(Tabij, k);
    energyMoments[k-1] = evaluateMp2(mp2Contribution.str(), Tabij);
    // Hartree--Fock correction (1st FT-MBPT)
    Tensor<> Tij(2, &Tabij.lens[2]);
    computeThermalHfAmplitudes(Tij, k);
    std::stringstream hfContribution;
    hfContribution << k << ".central moment of H^HF";
    // TODO: H_0 contribution: log(Z_0(beta)) = sum_1^Np log(1+exp(-eps_p*beta))
    energyMoments[k-1] += evaluateHf(hfContribution.str(), Tij);
  }
  std::vector<int64_t> indices(Dabij->wrld->rank == 0 ? n : 0);
  for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; };
  auto ctfEnergyMoments(new CTF::Tensor<>(1, std::vector<int>({int(n)}).data()) );
  ctfEnergyMoments->write(indices.size(), indices.data(), energyMoments.data());
  allocatedTensorArgument("ThermalMp2EnergyMoments", ctfEnergyMoments);
}

/**
 * \brief Computes the nth derivative of the thermal MP2 amplitudes
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * and returns the result in the tensor Tabij.
 * Expects Tabij to be zero on entry.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::computeThermalMp2Amplitudes(
  Tensor<> &Tabij, const int n
) {
  if (n == 0) {
    // all derivative degrees are zero
    addThermalMp2Amplitudes(Tabij, std::vector<int>(5));
  } else {
    LOG(1, "FT-MP2") << "computing MP2 contribution of d^"
      << n << " log(Z(beta)) / d(-beta)^" << n << " ..." << std::endl;
    // go through all possiblities to distribute n derivatives among 5 terms
    for (auto &derivedTerms: MultiCombinations(5, n)) {
      // convert to degree of derivative for each term
      std::vector<int> degrees(5);
      for (auto term: derivedTerms) {
        ++degrees[term];
      };
      addThermalMp2Amplitudes(Tabij, degrees);
    }
  }
}

/**
 * \brief Adds the contribution of the thermal MP2 amplitudes
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * to the tensor Tabij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 5 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addThermalMp2Amplitudes(
  Tensor<> &Tabij, const std::vector<int> &degrees
) {
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  // start with Vabij
  Tensor<> tabij(*Vabij);
  // apply the derivative of the 5 terms of the respective degrees:
  // Tabij *= f^a = 1/(1+exp(-eps_a*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[0])
    )
  ) (
    (*epsa)["a"], tabij["abij"]
  );
  // Tabij *= f^b
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[1])
    )
  ) (
    (*epsa)["b"], tabij["abij"]
  );
  // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[2])
    )
  ) (
    (*epsi)["i"], tabij["abij"]
  );
  // Tabij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[3])
    )
  ) (
    (*epsi)["j"], tabij["abij"]
  );
  // Tabij *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalMp2Propagation<>(beta, degrees[4])
    )
  ) (
    (*Dabij)["abij"], tabij["abij"]
  );
  // add this contribution
  Tabij["abij"] += tabij["abij"];
}

real ThermalMp2EnergyFromCoulombIntegrals::evaluateMp2(
  const std::string &contribution,
  Tensor<> &Tabij,
  real f
) {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  energy[""] = f * (+0.5) * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real direct( energy.get_val() );
  energy[""] = f *(-0.5) * spins * Tabij["abij"] * (*Vabij)["abji"];
  real exchange( energy.get_val() );

  LOG(0, "FT-MP2") << contribution << "=" << direct+exchange << std::endl;
  LOG(1, "FT-MP2") << contribution << "_d=" << direct << std::endl;
  LOG(1, "FT-MP2") << contribution << "_x=" << exchange << std::endl;

  return direct+exchange;
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ (-1) \int_0^\beta{\rm d}\tau_1\,f_i\, f_j \f$
 * and returns the result in the tensor Tij.
 * Expects Tij to be zero on entry.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::computeThermalHfAmplitudes(
  Tensor<> &Tij, const int n
) {
  if (n == 0) {
    // all derivative degrees are zero
    addThermalHfAmplitudes(Tij, std::vector<int>(3));
  } else {
    LOG(1, "FT-MP2") << "computing HF contribution of d^"
      << n << " log(Z(beta)) / d(-beta)^" << n << " ..." << std::endl;
    // go through all possiblities to distribute n derivatives among 3 terms
    for (auto &derivedTerms: MultiCombinations(3, n)) {
      // convert to degree of derivative for each term
      std::vector<int> degrees(3);
      for (auto term: derivedTerms) {
        ++degrees[term];
      };
      addThermalHfAmplitudes(Tij, degrees);
    }
  }
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ \int_0^\beta{\rm d}\tau\, f_i\, f_j \f$
 * to the tensor Tij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 3 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addThermalHfAmplitudes(
  Tensor<> &Tij, const std::vector<int> &degrees
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> tij(false, Tij);
  // start with 1

  // the first term is (-1) int_0^beta dtau = -beta
  switch (degrees[0]) {
  case 0:
    tij["ij"] = -beta;
    break;
  case 1:
    tij["ij"] = 1;
    break;
  default:
    // all higher derivatives of -beta w.r.t. (-beta) are zero
    // add nothing to Tabij and return
    return;
  }
  // apply the derivative of the other 2 terms of the respective degrees:
  // Tij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[1])
    )
  ) (
    (*epsi)["i"], tij["ij"]
  );
  // Tij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[2])
    )
  ) (
    (*epsi)["j"], tij["ij"]
  );
  // add this contribution
  Tij["ij"] += tij["ij"];
}

real ThermalMp2EnergyFromCoulombIntegrals::evaluateHf(
  const std::string &contribution,
  Tensor<> &Tij,
  real f
) {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  energy[""] = f * (+0.5) * spins * spins * Tij["ij"] * (*Vijkl)["ijij"];
  real direct( energy.get_val() );
  energy[""] = f * (-0.5) * spins * Tij["ij"] * (*Vijkl)["ijji"];
  real exchange( energy.get_val() );

  LOG(0, "FT-MP2") << contribution << "=" << direct+exchange << std::endl;
  LOG(1, "FT-MP2") << contribution << "_d=" << direct << std::endl;
  LOG(1, "FT-MP2") << contribution << "_x=" << exchange << std::endl;

  return direct+exchange;
}

