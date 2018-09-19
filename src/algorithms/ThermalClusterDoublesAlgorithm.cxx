#include <algorithms/ThermalClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <mixers/Mixer.hpp>
#include <util/LapackMatrix.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ThermalClusterDoublesAlgorithm::ThermalClusterDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalClusterDoublesAlgorithm::~ThermalClusterDoublesAlgorithm() {
}

void ThermalClusterDoublesAlgorithm::run() {
  beta = 1 / getRealArgument("Temperature");
  LOG(1, getCapitalizedAbbreviation()) << "beta=" << beta << std::endl;

  computeSqrtOccupancies();
  diagonalizeSinglesHamiltonian();

  real tda( getTammDancoffEnergy() );
  LOG(1, getCapitalizedAbbreviation()) << "TDA F=" << tda << std::endl;
  if (isArgumentGiven("ThermalTdaEnergy")) {
    setRealArgument("ThermalTdaEnergy", tda);
  }
  if (getIntegerArgument("linearized", 0)) return;


  // compute the other contributions perturbatively
  // get imaginary time and frequency grids on all nodes
  auto tn( getTensorArgument<real>("ImaginaryTimePoints") );
  taus.resize(tn->lens[0]);
  tn->read_all(taus.data());
  // truncate at beta
  while (taus.back() > beta) taus.pop_back();
  real q( getRealArgument("imaginaryTimeScale", 2.0) );
  // fill up with a geometric sequence until beta
  while (taus.back() < beta) taus.push_back(q*taus.back());
  taus.back() = beta;

  // get steady state amplitudes at tau points
//  computeSteadyStateAmplitudes(taus);
  if (getIntegerArgument("zeroTemperatureLimit", 0)) {
    real energy(getZeroTDrccd());
    std::stringstream energyName;
    energyName << "Thermal" << getAbbreviation() << "Energy";
    setRealArgument(energyName.str(), energy);
    return;
  }

//  iterateAmplitudeFunctions();
  iterateAmplitudeSamples();

  // diagonalizeDoublesAmplitudes();
}

void ThermalClusterDoublesAlgorithm::dryRun() {
}

std::string ThermalClusterDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(),
    abbreviation.begin(), ::toupper
  );
  return abbreviation;
}

std::string ThermalClusterDoublesAlgorithm::getAmplitudeIndices(Tensor<> &T) {
  char indices[T.order+1];
  const int excitationLevel(T.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    indices[i] = static_cast<char>('a'+i);
    indices[i+excitationLevel] = static_cast<char>('i'+i);
  }
  indices[T.order] = 0;
  return indices;
}

void ThermalClusterDoublesAlgorithm::iterateAmplitudeSamples() {
  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );

  std::vector<real> amplitudeEnergyDensities(taus.size());
  real energy;
  real accuracy(getRealArgument("accuracy", 1e-8));

  Tensor<real> T0FG(false, *VdFG);
  std::vector<PTR(Tensor<real>)> tensors({NEW(Tensor<real>, false, *VdFG)});
  std::vector<std::string> indices({"FG"});
  PTR(const FockVector<real>) amplitudes(
    NEW(FockVector<real>,
      tensors.begin(), tensors.end(),
      indices.begin(), indices.end()
    )
  );
  real DTau, tau0(0.0);
  // solve amplitudes on imaginary time grid one-by-one
  real direct(0.0), exchange(0.0);
  real d(0.0), x(0.0);
  for (size_t n(0); n < taus.size(); ++n) {
    real tau1(taus[n]);
    DTau = tau1-tau0;
    // add energy contribution from previous interval
    direct += 0.5*d*DTau; exchange += 0.5*x*DTau;
    LOG(1, getCapitalizedAbbreviation())
      << "solving amplitudes at tau_" << (n+1) << "=" << taus[n]
      << std::endl;
    real lastEnergy(0.0);
    // create a mixer, by default use the linear one
    std::string mixerName(getTextArgument("mixer", "LinearMixer"));
    Mixer<real> *mixer( MixerFactory<real>::create(mixerName, this) );
    if (!mixer) {
      std::stringstream stringStream;
      stringStream << "Mixer not implemented: " << mixerName;
      throw new EXCEPTION(stringStream.str());
    }
    for (int i(0); i < I; ++i) {
      LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
      auto estimatedAmplitudes( NEW(FockVector<real>, *amplitudes) );
      auto T1FG(amplitudes->get(0));
      auto S1FG(estimatedAmplitudes->get(0));
      (*S1FG)["FG"] = T0FG["FG"];
      // propagate amplitudes T0 from tau_n-1 to tau_n with effective H0
      ImaginaryTimePropagation imaginaryTimePropagation(DTau);
      propagateAmplitudes(*S1FG, imaginaryTimePropagation);
      // apply hamiltonian between tau_n-1 and tau_n
      applyHamiltonian(T0FG, *T1FG, DTau, *S1FG);
      // compute amplitudes change and tell mixer
      auto amplitudesChange( NEW(FockVector<real>, *estimatedAmplitudes) );
      *amplitudesChange -= *amplitudes;
      mixer->append(estimatedAmplitudes, amplitudesChange);
      // get mixer's best guess for amplitudes
      amplitudes = mixer->get();
      T1FG = amplitudes->get(0);

      d = 0.0; x = 0.0;
      computeEnergyContribution(*T1FG, 2.0, d, x);
      LOG(2, getCapitalizedAbbreviation()) << "F_d=" <<
        (direct+0.5*d*DTau)/tau1 << std::endl;
      LOG(2, getCapitalizedAbbreviation()) << "F_x=" <<
        (exchange+0.5*x*DTau)/tau1 << std::endl;
      // write energy densities
      amplitudeEnergyDensities[n] = d+x;
      LOG(2, getCapitalizedAbbreviation())
        << "T.V=" << amplitudeEnergyDensities[n] << std::endl;
      energy = direct + exchange + 0.5*(d+x)*DTau;
      LOG(1, getCapitalizedAbbreviation()) << "F=" << energy/tau1 << std::endl;
      if (std::abs(1-lastEnergy/energy) < accuracy) break;
      lastEnergy = energy;
    }
    // add energy contribution from this interval
    direct += 0.5*d*DTau; exchange += 0.5*x*DTau;
    // go to next interval
    tau0 = tau1;
    T0FG["FG"] = (*amplitudes->get(0))["FG"];
  }

  LOG(1, getCapitalizedAbbreviation()) << "F_d=" << direct/beta << std::endl;
  LOG(1, getCapitalizedAbbreviation()) << "F_x=" << exchange/beta << std::endl;

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), energy/beta);

  if (isArgumentGiven("Amplitudes")) {
    auto newTFG(new CTF::Tensor<real>(*amplitudes->get(0)));
    allocatedTensorArgument<real>("Amplitudes", newTFG);
  }
  // write amplitude norms
  std::vector<int64_t> normIndices(T0FG.wrld->rank == 0 ? taus.size() : 0);
  for (size_t i(0); i < normIndices.size(); ++i) { normIndices[i] = i; }
  if (isArgumentGiven("amplitudeTimes")) {
    int n(taus.size());
    auto taun(new Tensor<>(1, &n, T0FG.sym, *T0FG.wrld, "amplitudeTimes"));
    taun->write(
      normIndices.size(), normIndices.data(), taus.data()
    );
    allocatedTensorArgument<real>("amplitudeTimes", taun);
  }
  if (isArgumentGiven("amplitudeNorms")) {
    int n(taus.size());
    auto normn(new Tensor<>(1, &n, T0FG.sym, *T0FG.wrld, "amplitudeNorms"));
    normn->write(
      normIndices.size(), normIndices.data(), amplitudeEnergyDensities.data()
    );
    allocatedTensorArgument<real>("amplitudeNorms", normn);
  }
}

void ThermalClusterDoublesAlgorithm::iterateAmplitudeFunctions() {
  auto Vabij( getTensorArgument<real>("ThermalPPHHCoulombIntegrals") );

  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  // allocate doubles amplitudes on imaginary time grid
  // the source amplitudes are stored in orbital space
  // with fully closed contraction weights, suitable for interpolation
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }
  // the target amplitudes
  std::vector<PTR(Tensor<real>)> Sabijn(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Sabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }
  std::vector<real> energies(taus.size());

  // TODO: test then accept/dismiss renormalizations, for now disabled
  int R( getIntegerArgument("renormalizations", 0) );

  real energy;
  real accuracy(getRealArgument("accuracy", 1e-7));
  CTF::Tensor<real> T0abij(false, *Vabij);
  CTF::Tensor<real> S1abij(false, *Vabij);
  for (int r(-R); r <= 0; ++r) {
    real scale( std::pow(taus.back()/taus.front(),r) );
    LOG(0, getCapitalizedAbbreviation()) << "renormalization level: " << r << std::endl;
    {
      size_t N(taus.size()-1);
//    for (size_t N(0); N < taus.size(); ++N) {
      real lastEnergy(0);
      for (int i(0); i < I; ++i) {
        LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
        real tau0(0.0);
        T0abij["abij"] = 0.0;
        S1abij["abij"] = 0.0;
        real direct(0.0), exchange(0.0);
        for (size_t n(0); n <= N; ++n) {
          real tau1(scale*taus[n]);
          real DTau(tau1-tau0);
          // energy contribution from previously convolved amplitudes ST^I(tau_n-1)
          computeEnergyContribution(S1abij, DTau, direct, exchange);

          // get T0abij at tau_n-1 and write previously convolved amplitudes
          // note T(0) is implicitly 0
          if (n > 0) {
            T0abij["abij"] = (*Tabijn[n-1])["abij"];
            (*Sabijn[n-1])["abij"] = S1abij["abij"];
          }

          LOG(1, getCapitalizedAbbreviation())
            << "convolving amplitudes at tau_" << (n+1) << "=" << tau1
            << std::endl;

          // propagate previously convolved amplitudes ST^I(tau_n-1) to this tau_n
          ImaginaryTimePropagation imaginaryTimePropagation(DTau);
          propagateAmplitudes(S1abij, imaginaryTimePropagation);

          // apply hamiltonian between tau_n-1 and tau_n
          applyHamiltonian(T0abij, *Tabijn[n], DTau, S1abij);

          // energy contribution from convolved amplitudes ST^I(tau_n)
          computeEnergyContribution(S1abij, DTau, direct, exchange);
          LOG(2, getCapitalizedAbbreviation())
            << "F_d=" << direct/tau1 << std::endl;
          LOG(2, getCapitalizedAbbreviation())
            << "F_x=" << exchange/tau1 << std::endl;
          // write norm
          real a(S1abij.norm2());
          LOG(2, getCapitalizedAbbreviation()) << "|T|=" << a << std::endl;
          tau0 = tau1;
        }
        (*Sabijn[N])["abij"] = S1abij["abij"];
        energies[N] = energy = direct + exchange;
        LOG(1, getCapitalizedAbbreviation()) << "F=" << energy/tau0 << std::endl;
        if (std::abs(1-lastEnergy/energy) < accuracy) break;
        lastEnergy = energy;
        // mix with previous step
        real mixingRatio( getRealArgument("mixingRatio", 1.0) );
        for (size_t n(0); n <= N; ++n) {
          (*Tabijn[n])["abij"] *= (1-mixingRatio);
          // mix
          (*Tabijn[n])["abij"] += mixingRatio * (*Sabijn[n])["abij"];
        }
      }
    }
  }

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), energy/beta);

  if (isArgumentGiven("Amplitudes")) {
    auto newTabij(new CTF::Tensor<real>(S1abij));
    allocatedTensorArgument<real>("Amplitudes", newTabij);
  }
}

cc4s::real ThermalClusterDoublesAlgorithm::getTammDancoffEnergy() {
  real spins( getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0 );

  // compute Tamm-Dancoff Approximation (TDA)
  Scalar<> e;
  int NF(lambdaF->lens[0]);

  // propagate doubles
  lambdaFG = NEW(Tensor<real>, false, *VdFG);
  (*lambdaFG)["FG"] =  (*lambdaF)["F"];
  (*lambdaFG)["FG"] += (*lambdaF)["G"];
  Tensor<real> TFG(*VdFG);
  SecondOrderIntegral secondOrderIntegral(beta);
  Transform<real, real>(std::function<void(real,real&)>(secondOrderIntegral))(
    (*lambdaFG)["FG"], TFG["FG"]
  );
  // project out null-space from Fock-space-doubling
  Transform<real, real> projectOut(
    std::function<void(real, real &)>(
      [](const real lambda, real &t) {
        if (lambda == 0.0) t = 0;
      }
    )
  );
  projectOut( (*lambdaF)["F"], TFG["FG"] );
  projectOut( (*lambdaF)["G"], TFG["FG"] );
  e[""] = -0.5 * spins*spins * TFG["FG"] * (*VdFG)["FG"];
  real tdad(e.get_val());
  LOG(1, getCapitalizedAbbreviation()) << "TDA F_d=" << tdad << std::endl;

  e[""] = +0.5 * spins * TFG["FG"] * (*VxFG)["FG"];
  real tdax(e.get_val());
  LOG(1, getCapitalizedAbbreviation()) << "TDA F_x=" << tdax << std::endl;

  if (isArgumentGiven("SinglesHamiltonianWeights")) {
    auto singlesWeights(new Tensor<real>(1, std::vector<int>({NF}).data()) );
    (*singlesWeights)["F"] = 0.5 * spins*spins * TFG["FF"] * (*VdFG)["FF"];
    allocatedTensorArgument<real>("SinglesHamiltonianWeights", singlesWeights);
  }

  return tdad+tdax;
}

void ThermalClusterDoublesAlgorithm::computeEnergyContribution(
  Tensor<real> &SFG, const real DTau,
  real &direct, real &exchange
) {
//  auto Vabij( getTensorArgument<>("ThermalPPHHCoulombIntegrals") );
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<real> energy;

  // direct term
  energy[""] = (+0.5) * spins * spins * DTau/2 * SFG["FG"] * (*VdFG)["FG"];
//  energy[""] = 0.5 * spins * spins * DTau/2 * Sabij["abij"] * (*Vabij)["abij"];
//    * (*ga)["a"] * (*ga)["b"] * (*gi)["i"] * (*gi)["j"];
  direct += energy.get_val();
  // exchange term
  energy[""] = (-0.5) * spins * DTau/2 * SFG["FG"] * (*VxFG)["FG"];
//  energy[""] = (-0.5) * spins * DTau/2 * Sabij["abij"] * (*Vabij)["abji"];
//    * (*ga)["a"] * (*ga)["b"] * (*gi)["i"] * (*gi)["j"];
  exchange += energy.get_val();
}

void ThermalClusterDoublesAlgorithm::computeSqrtOccupancies() {
  // compute g^a and g_i = sqrt(occupancies)
  auto Ni( getTensorArgument<>("ThermalHoleOccupancies") );
  auto Na( getTensorArgument<>("ThermalParticleOccupancies") );
  gi = NEW(Tensor<real>, *Ni);
  Transform<real>(
    std::function<void(real &)>([](real &f) { f = std::sqrt(f); } )
  ) (
    (*gi)["i"]
  );
  ga = NEW(Tensor<real>, *Na);
  Transform<real>(
    std::function<void(real &)>([](real &f) { f = std::sqrt(f); } )
  ) (
    (*ga)["a"]
  );
}

void ThermalClusterDoublesAlgorithm::diagonalizeSinglesHamiltonian() {
  auto epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  auto Vabij(getTensorArgument<>("ThermalPPHHCoulombIntegrals"));
  auto Vbija(getTensorArgument<>("ThermalPHHPCoulombIntegrals"));
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // build to Hbjai
  int Nv(epsa->lens[0]); int No(epsi->lens[0]);
  int lens[] = { Nv,No, Nv,No };
  auto Hbjai(NEW(Tensor<>, 4, lens, Vbija->sym, *Vbija->wrld, "Hbjai"));
  // bubble from H_1
  (*Hbjai)["bjai"] = spins * (*Vbija)["bija"];

  // half-close all contractions on inserted interactions
  (*Hbjai)["bjai"] *= (*ga)["b"];
  (*Hbjai)["bjai"] *= (*gi)["j"];
  (*Hbjai)["bjai"] *= (*ga)["a"];
  (*Hbjai)["bjai"] *= (*gi)["i"];

  // unperturbed propatation is diagonal
  (*Hbjai)["bjbj"] += (*epsa)["b"];
  (*Hbjai)["bjbj"] -= (*epsi)["j"];

  LOG(1, getCapitalizedAbbreviation())
    << "diagonalizing singles part of Hamiltonian..." << std::endl;
  // H(bj)(ai) = U.S.U^T, seen as a matrix with compound indices
  BlacsWorld world(Hbjai->wrld->rank, Hbjai->wrld->np);
  int NvNo(Nv*No);
  int scaHLens[2] = { NvNo, NvNo };
  auto scaHbjai(NEW(ScaLapackMatrix<>, *Hbjai, scaHLens, &world));
  // release unneeded resources early
  Hbjai = nullptr;

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(NEW(ScaLapackMatrix<>, *scaHbjai));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaHbjai, scaU);
  std::vector<real> lambdas(NvNo);
  eigenSystem.solve(lambdas.data());
  scaHbjai = nullptr;

  // determine dim(ker(effective Hamiltonian)) from Fock-space doubling
  std::vector<real> occ(No), eps(No);
  epsi->read_all(eps.data());
  getTensorArgument<real>("ThermalHoleOccupancies")->read_all(occ.data());
  size_t dimKerH0(0), degenerateStart(0);
  real degenerateEps(eps[0]);
  real degeneracyThreshold(getRealArgument("degeneracyThreshold", 1e-6));
  for (size_t i(0); i < eps.size(); ++i) {
    if (eps[i] - degenerateEps > degeneracyThreshold) {
      if (occ[degenerateStart] < 1.0) {
        size_t n(i-degenerateStart);
        dimKerH0 += n*(n-1)/2;
        if (n > 1) {
          LOG(2, getCapitalizedAbbreviation()) <<
            n << " overlapping degenerate states at epsilon=" << degenerateEps << std::endl;
        }
      }
      degenerateEps = eps[i];
      degenerateStart = i;
    }
  }
  if (occ[degenerateStart] < 1.0) {
    size_t n(eps.size()-degenerateStart);
    dimKerH0 += n*(n-1)/2;
    if (n > 1) {
      LOG(2, getCapitalizedAbbreviation()) <<
        n << " overlapping degenerate states at epsilon=" << degenerateEps << std::endl;
    }
  }
//  dimKerH0 += 13;
  LOG(1, getCapitalizedAbbreviation())
    << "dim(ker(H_singles))=" << dimKerH0 << std::endl;

  // get lowest |eigenvalues| from spectrum
  std::list<std::pair<size_t,real>> lowestNegativeEigenvalues;
  size_t F(0);
  while (F < lambdas.size()) {
    if (lambdas[F] < 0.0) {
      lowestNegativeEigenvalues.push_front(std::make_pair(F,-lambdas[F]));
    } else {
      break;
    }
    ++F;
  }
  real modeThreshold(1.0);
  size_t nullModesCount(0);
  size_t nullSpaceStart(0), nullSpaceEnd(0);
  while (nullModesCount < dimKerH0) {
    if (
      lowestNegativeEigenvalues.size() > 0 &&
      lowestNegativeEigenvalues.front().second < lambdas[F]
    ) {
      nullSpaceStart = lowestNegativeEigenvalues.front().first;
      modeThreshold = lowestNegativeEigenvalues.front().second;
      lowestNegativeEigenvalues.pop_front();
    } else {
      nullSpaceEnd = F+1;
      modeThreshold = lambdas[F];
      ++F;
    }
    ++nullModesCount;
  }
  real minMode(0.0);
  if (
    lowestNegativeEigenvalues.size() > 0 &&
    lowestNegativeEigenvalues.front().second < lambdas[F]
  ) {
    minMode = lowestNegativeEigenvalues.front().second;
  } else {
    minMode = lambdas[F];
  }
  LOG(1, getCapitalizedAbbreviation()) <<
    "nullspace modes " << nullSpaceStart << "<=F<" << nullSpaceEnd << std::endl;
  if (nullSpaceStart < nullSpaceEnd) {
    LOG(1, getCapitalizedAbbreviation()) <<
      "largest null mode " << modeThreshold << std::endl;
    LOG(1, getCapitalizedAbbreviation()) <<
      "smallest non-null mode " << minMode << std::endl;
    for (size_t F(nullSpaceStart); F < nullSpaceEnd; ++F) lambdas[F] = 0.0;
  }

  // write Lambda back to CTF
  std::vector<int64_t> lambdaIndices(epsi->wrld->rank == 0 ? NvNo : 0);
  for (size_t i(0); i < lambdaIndices.size(); ++i) { lambdaIndices[i] = i; }
  lambdaF = new Tensor<>(1, &NvNo, Vbija->sym, *Vbija->wrld, "Lambda");
  lambdaF->write(lambdaIndices.size(), lambdaIndices.data(), lambdas.data());

  allocatedTensorArgument<>("SinglesHamiltonianEigenvalues", lambdaF);

  // write unitary matrix U(ai)(F) back to CTF as tensor UaiF
  int ULens[3] = { Nv, No, NvNo };
  UaiF = NEW(Tensor<>, 3, ULens, Vbija->sym, *Vbija->wrld, "UaiF");
  scaU->write(*UaiF);
  scaU = nullptr;
  // chop near-zero values to zero of unitary transform
  Transform<real> chop(
    std::function<void(real &)>(
      [](real &t) {
        if (std::abs(t) < sqrt(std::numeric_limits<real>::epsilon())) t=0.0;
      }
    )
  );
  chop((*UaiF)["aiF"]);
  if (
    nullSpaceStart < nullSpaceEnd &&
    isArgumentGiven("singlesHamiltonianNullspace")
  ) {
    int UZerosStart[] = {0, 0, static_cast<int>(nullSpaceStart)};
    int UZerosEnd[] = {Nv, No, static_cast<int>(nullSpaceEnd)};
    allocatedTensorArgument<real>(
      "singlesHamiltonianNullspace",
      new Tensor<real>(UaiF->slice(UZerosStart, UZerosEnd))
    );
  }

  // Coulomb interaction in singles-mode space, half-closed
  int NF(lambdaF->lens[0]);
  // two particle/hole pairs F&G
  VdFG = NEW(Tensor<real>, 2, std::vector<int>({NF,NF}).data());
  (*VdFG)["FG"] = (*UaiF)["ckF"] * (*UaiF)["dlG"] *
    (*ga)["c"] * (*ga)["d"] * (*gi)["k"] * (*gi)["l"] *
    (*Vabij)["cdkl"];
  chop((*VdFG)["FG"]);
  // exchange interaction
  VxFG = NEW(Tensor<real>, false, *VdFG);
  (*VxFG)["FG"] = (*UaiF)["ckF"] * (*UaiF)["dlG"] *
    (*ga)["c"] * (*ga)["d"] * (*gi)["k"] * (*gi)["l"] *
    (*Vabij)["cdlk"];
  chop((*VxFG)["FG"]);
  // determine Coulomb coupling to nullspapce
  if (nullSpaceStart < nullSpaceEnd) {
    int VZerosStart[] = {0, static_cast<int>(nullSpaceStart)};
    int VZerosEnd[] = {NvNo, static_cast<int>(nullSpaceEnd)};
    Tensor<real> coupling(VdFG->slice(VZerosStart, VZerosEnd));
    Vector<real> norms(VZerosEnd[1]-VZerosStart[1]);
    norms["F"] = coupling["GF"] * coupling["GF"];
    real couplingNorm(std::sqrt(norms.norm_infty()));
    LOG(1, getCapitalizedAbbreviation()) <<
      "max(|Coulomb coupling to nullspace|)=" << couplingNorm << std::endl;
    if (isArgumentGiven("couplingToSinglesNullspace")) {
      allocatedTensorArgument<real>(
        "couplingToSinglesNullspace", new Tensor<real>(coupling)
      );
    }
  }
  // truncate coupling to nullspace
  Transform<real, real> projectOut(
    std::function<void(real, real &)>(
      [](const real lambda, real &t) {
        if (lambda == 0.0) t = 0;
      }
    )
  );
  projectOut( (*lambdaF)["F"], (*VdFG)["FG"] );
  projectOut( (*lambdaF)["G"], (*VdFG)["FG"] );
}

real ThermalClusterDoublesAlgorithm::getZeroTDrccd() {
  real spins( getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0 );
  real levelShift( getRealArgument("levelShift", 0.0) );
  // level shifted division for left hand side
  class LevelShiftedDivision {
  public:
    LevelShiftedDivision(real shift_): shift(shift_) { }
    void operator()(real lambda, real &s) {
      if (lambda == 0.0) {
        s = 0;
      } else {
        s = -s / (lambda + shift);
      }
    }
  protected:
    real shift;
  } levelShiftedDivision(levelShift);

  // create a mixer, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  Mixer<real> *mixer( MixerFactory<real>::create(mixerName, this) );
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  std::vector<PTR(Tensor<real>)> tensors({NEW(Tensor<real>, false, *VdFG)});
  std::vector<std::string> indices({"FG"});
  PTR(const FockVector<real>) amplitudes(
    NEW(FockVector<real>,
      tensors.begin(), tensors.end(),
      indices.begin(), indices.end()
    )
  );
  Scalar<> e;
  real energy(0), lastEnergy(0);
  real accuracy(getRealArgument("accuracy", 1e-7));
  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  for (int i(0); i < I; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    real direct(0.0), exchange(0.0);
    LOG(1, getCapitalizedAbbreviation())
      << "solving for drCCD steady state amplitudes" << std::endl;
    auto estimatedAmplitudes( NEW(FockVector<real>, *amplitudes) );
    auto TFG(amplitudes->get(0));
    auto SFG(estimatedAmplitudes->get(0));
    Tensor<real> WFG(false, *VdFG);
    WFG["FG"]  = (+1.0) * spins * spins * (*VdFG)["FG"];
    if (getIntegerArgument("adjacentPairsExchange", 0)) {
      WFG["FG"] += (-1.0) * spins * (*VxFG)["FG"];
    } 
    // constant term
    (*SFG)["FG"]  = (*VdFG)["FG"];
    // quadratic term
    (*SFG)["FG"] += (*TFG)["FH"] * WFG["HI"] * (*TFG)["IG"];
    // apply level shifting on right hand side
    (*SFG)["FG"] += (-levelShift) * (*TFG)["FG"];
    Transform<real, real> projectOut(
      std::function<void(real, real &)>(
        [](const real lambda, real &t) {
          if (lambda == 0.0) t = 0;
        }
      )
    );
    projectOut( (*lambdaF)["F"], (*SFG)["FG"] );
    projectOut( (*lambdaF)["G"], (*SFG)["FG"] );
    // divide by -(Delta+shift) to get new estimate for T
    Transform<real, real>(
      std::function<void(real, real &)>(levelShiftedDivision)
    ) (
      (*lambdaFG)["FG"], (*SFG)["FG"]
    );
    // compute amplitudes change and tell mixer
    auto amplitudesChange( NEW(FockVector<real>, *estimatedAmplitudes) );
    *amplitudesChange -= *amplitudes;
    mixer->append(estimatedAmplitudes, amplitudesChange);
    // get mixer's best guess for amplitudes
    amplitudes = mixer->get();
    TFG = amplitudes->get(0);
    // write norm
    real l2((*TFG).norm2()), linf((*TFG).norm_infty());
    LOG(2, getCapitalizedAbbreviation())
      << "|T|=" << l2 << ", max(T)=" << linf << std::endl;
    e[""] = +0.5 * spins*spins * (*TFG)["FG"] * (*VdFG)["FG"];
    direct = e.get_val();
    LOG(1, getCapitalizedAbbreviation()) << "T->0 F_d=" << direct << std::endl;
    e[""] = -0.5 * spins * (*TFG)["FG"] * (*VxFG)["FG"];
    exchange = e.get_val();
    LOG(1, getCapitalizedAbbreviation()) << "T->0 F_x=" << exchange << std::endl;
    energy = direct + exchange;
    LOG(1, getCapitalizedAbbreviation()) << "T->0 F=" << energy << std::endl;
    if (std::abs(1-lastEnergy/energy) < accuracy) break;
    lastEnergy = energy;
  }
  return energy;
}

void ThermalClusterDoublesAlgorithm::computeSteadyStateAmplitudes(
  const std::vector<real> &taus
) {
  auto Vabij( getTensorArgument<>("ThermalPPHHCoulombIntegrals") );
  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  real energy;
  real accuracy(getRealArgument("accuracy", 1e-7));
  CTF::Tensor<real> T0abij(false, *Vabij);
  CTF::Tensor<real> S1abij(false, *Vabij);
  T0abij["abij"] = 0.0;
  real lastEnergy(0);
  for (size_t n(0); n < taus.size(); ++n) {
    LOG(0, getCapitalizedAbbreviation()) << "renormalization tau_" << (n+1)
      << "=" << taus[n] << std::endl;
    for (int i(0); i < I; ++i) {
      LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
      real direct(0.0), exchange(0.0);
      LOG(1, getCapitalizedAbbreviation())
        << "convolving steady state amplitudes" << std::endl;
      // apply hamiltonian between 0 and tau_n
      S1abij["abij"] = 0.0;
      applyHamiltonian(T0abij, T0abij, taus[n], S1abij);
      // energy contribution from convolved amplitudes ST^I(tau_n)
      computeEnergyContribution(S1abij, 2.0, direct, exchange);
      LOG(2, getCapitalizedAbbreviation()) << "F_d=" << direct << std::endl;
      LOG(2, getCapitalizedAbbreviation()) << "F_x=" << exchange << std::endl;
      // write norm
      real a(S1abij.norm2());
      LOG(2, getCapitalizedAbbreviation()) << "|T|=" << a << std::endl;
      energy = direct + exchange;
      LOG(1, getCapitalizedAbbreviation()) << "F=" << energy << std::endl;
      if (std::abs(1-lastEnergy/energy) < accuracy) break;
      lastEnergy = energy;
      real mixingRatio( getRealArgument("mixingRatio", 1.0) );
      T0abij["abij"] *= (1-mixingRatio);
      T0abij["abij"] += mixingRatio * S1abij["abij"];
    }
  }  
}

void ThermalClusterDoublesAlgorithm::propagateAmplitudes(
  Tensor<real> &SFG,
  const std::function<void(real, real &)> &propagator
) {
/*
  Transform<real> chop(
    std::function<void(real &)>(
      [](real &t) {
        if (std::abs(t) < sqrt(std::numeric_limits<real>::epsilon())) t = 0.0;
      }
    )
  );
  chop( Sabij["abij"] );
  Transform<real, real> divBy(
    std::function<void(real, real &)>([](const real g, real &t){ t /= g; })
  );
  // half-open contraction weights for propagation
  divBy( (*ga)["a"], Sabij["abij"] );
  divBy( (*ga)["b"], Sabij["abij"] );
  divBy( (*gi)["i"], Sabij["abij"] );
  divBy( (*gi)["j"], Sabij["abij"] );
  // transform to singles-mode space
  Tensor<real> SFG(false, *VdFG);
  SFG["FG"] = (*UaiF)["bjG"] * (*UaiF)["aiF"] * Sabij["abij"];
*/
  // propagate
  Transform<real, real>(std::function<void(real, real &)>(propagator))(
    (*lambdaFG)["FG"], SFG["FG"]
  );
  // project-out Fock space doubling DOF in effective hamiltonian
  Transform<real, real> projectOut(
    std::function<void(real, real &)>(
      [](const real lambda, real &t) {
        if (lambda == 0.0) t = 0;
      }
    )
  );
  projectOut( (*lambdaF)["F"], SFG["FG"] );
  projectOut( (*lambdaF)["G"], SFG["FG"] );
/*
  // transform back to orbital space and close contraction weights
  Sabij["abij"] = (*ga)["a"] * (*ga)["b"] * (*gi)["i"] * (*gi)["j"] *
    (*UaiF)["aiF"] * (*UaiF)["bjG"] * SFG["FG"];
*/
}

void ThermalClusterDoublesAlgorithm::diagonalizeDoublesAmplitudes() {
  LOG(1, getCapitalizedAbbreviation())
    << "diagonalizing doubles amplitudes at tau_1..." << std::endl;
  BlacsWorld world(Tabijn.back()->wrld->rank, Tabijn.back()->wrld->np);
  auto epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  int NvNo(epsa->lens[0]*epsi->lens[0]);
  int scaTLens[2] = { NvNo, NvNo };
  auto scaTFG(NEW(ScaLapackMatrix<>, *Tabijn.back(), scaTLens, &world));

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(NEW(ScaLapackMatrix<>, *scaTFG));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTFG, scaU);
  std::vector<real> lambdas(Tabijn.back()->lens[0]);
  eigenSystem.solve(lambdas.data());
  scaTFG = nullptr;

  if (isArgumentGiven("DoublesAmplitudesEigenvalues")) {
    // write Lambda back to CTF
    std::vector<int64_t> lambdaIndices(
      UaiF->wrld->rank == 0 ? lambdas.size() : 0
    );
    for (size_t i(0); i < lambdaIndices.size(); ++i) { lambdaIndices[i] = i; }
    auto lambdaL(
      new Tensor<>(
        1, scaTLens, Tabijn.front()->sym, *Tabijn.front()->wrld, "L"
      )
    );
    lambdaL->write(lambdaIndices.size(), lambdaIndices.data(), lambdas.data());
    allocatedTensorArgument<>("DoublesAmplitudesEigenvalues", lambdaF);
  }
}

