#include <limits>

#include <algorithms/CoulombVertexDecomposition.hpp>
#include <math/CanonicalPolyadicDecomposition.hpp>
#include <math/IterativePseudoInverse.hpp>
#include <math/RandomTensor.hpp>
#include <math/MathFunctions.hpp>
#include <mixers/Mixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

DEFSPEC(
    CoulombVertexDecomposition,
    SPEC_IN(
        {"delta", SPEC_VALUE_DEF("TODO: DOC", double, DEFAULT_DELTA)},
        {"rankFactor",
         SPEC_VALUE_DEF("TODO: DOC", double, DEFAULT_RANK_FACTOR)},
        {"swampingThreshold",
         SPEC_VALUE_DEF("TODO: DOC", double, DEFAULT_SWAMPING_THRESHOLD)},
        {"fitCoulombFactors", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"fitFactorOrbitals", SPEC_VALUE_DEF("TODO: DOC", int64_t, 1)},
        {"maxIterations",
         SPEC_VALUE_DEF("TODO: DOC", int64_t, DEFAULT_MAX_ITERATIONS)},
        {"maxSubIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 8)},
        {"minSubIterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 2)},
        {"rankSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, DEFAULT_RANK_SIZE)},
        {"realFactorOrbitals",
         SPEC_VALUE_DEF("TODO: DOC", int64_t, DEFAULT_REAL_FACTOR_ORBITALS)},
        {"writeSubIterations",
         SPEC_VALUE_DEF("TODO: DOC", int64_t, DEFAULT_WRITE_SUB_ITERATIONS)},
        {"ansatz", SPEC_VALUE_DEF("TODO: DOC", std::string, HERMITIAN)},
        {"mixer", SPEC_VALUE_DEF("TODO: DOC", std::string, "LinearMixer")},
        {"CoulombVertex", SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},
        {"StartingCoulombFactors",
         SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)},
        {"StartingFactorOrbitals",
         SPEC_VARIN("TODO: DOC", Tensor<sisi4s::complex> *)}),
    SPEC_OUT(
        {"ComposedCoulombVertex",
         SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)},
        {"CoulombFactors", SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)},
        {"FactorOrbitals", SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)},
        {"OutgoingFactorOrbitals",
         SPEC_VAROUT("TODO: DOC", Tensor<sisi4s::complex> *)}));

IMPLEMENT_ALGORITHM(CoulombVertexDecomposition) {
  GammaGqr = in.get<Tensor<sisi4s::complex> *>("CoulombVertex");
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);

  // calculate decomposition rank
  rank = in.get<int64_t>("rankSize", DEFAULT_RANK_SIZE);
  // if rank is not given use rank factors (if they are not given use
  // rankFactors=3.0)
  if (rank == -1) {
    double rankFactor(in.get<double>("rankFactor", DEFAULT_RANK_FACTOR));
    rank = NG * rankFactor;
  }

  realFactorOrbitals =
      in.get<int64_t>("realFactorOrbitals", DEFAULT_REAL_FACTOR_ORBITALS);
  normalizedFactorOrbitals =
      in.get<int64_t>("normalizedFactorOrbitals",
                      DEFAULT_NORMALIZED_FACTOR_ORBITALS);
  LOG(0, "RALS") << "Tensor rank decomposition with rank NR=" << rank
                 << ", realFactorOrbitals=" << realFactorOrbitals
                 << ", normalizedFactorOrbitals=" << normalizedFactorOrbitals
                 << std::endl;
  LOG(1, "RALS") << "Decomposing Coulomb vertex " << GammaGqr->get_name()
                 << " with NG=" << NG << ", Np=" << Np << std::endl;

  writeSubIterations =
      in.get<int64_t>("writeSubIterations", DEFAULT_WRITE_SUB_ITERATIONS);

  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);

  // allocate factor tensors
  if (isArgumentGiven("StartingFactorOrbitals")) {
    Tensor<sisi4s::complex> *PirRTensor(
        in.get<Tensor<sisi4s::complex> *>("StartingFactorOrbitals"));
    PirRTensor->set_name("StartingPirR");
    if (PirRTensor->order != 2)
      throw new EXCEPTION("Matrix expected as argument StartingPirR");
    LOG(1, "RALS") << "Initial PirR=" << PirRTensor->get_name() << std::endl;
    PirR = reinterpret_cast<Matrix<complex> *>(PirRTensor);
  } else {
    PirR = new Matrix<complex>(Np,
                               int(rank),
                               NS,
                               *GammaGqr->wrld,
                               "PirR",
                               GammaGqr->profile);
    LOG(1, "RALS") << "Initial PirR=RandomTensor" << std::endl;
    setRandomTensor(*PirR, normalDistribution, random);
    realizePi(*PirR);
    normalizePi(*PirR);
  }

  if (isArgumentGiven("StartingCoulombFactors")) {
    Tensor<sisi4s::complex> *LambdaGRTensor(
        in.get<Tensor<sisi4s::complex> *>("StartingCoulombFactors"));
    LambdaGRTensor->set_name("StartingLambdaGR");
    if (LambdaGRTensor->order != 2)
      throw new EXCEPTION("Matrix expected as argument StartingLambdaGR");
    LOG(1, "RALS") << "Initial LambdaGR=" << LambdaGRTensor->get_name()
                   << std::endl;
    LambdaGR = reinterpret_cast<Matrix<complex> *>(LambdaGRTensor);
  } else {
    LambdaGR = new Matrix<complex>(NG,
                                   int(rank),
                                   NS,
                                   *GammaGqr->wrld,
                                   "LambdaGR",
                                   GammaGqr->profile);
    LOG(1, "RALS") << "Initial LambdaGR=RandomTensor" << std::endl;
    setRandomTensor(*LambdaGR, normalDistribution, random);
  }

  PiqR = new Matrix<complex>(Np,
                             int(rank),
                             NS,
                             *GammaGqr->wrld,
                             "PiqR",
                             GammaGqr->profile);
  std::string ansatz(in.get<std::string>("ansatz", HERMITIAN));
  if (ansatz != HERMITIAN && ansatz != SYMMETRIC && ansatz != PSEUDO_INVERSE) {
    std::stringstream stringStream;
    stringStream << "Unknown decomposition ansatz \"" << ansatz << "\"";
    throw new EXCEPTION(stringStream.str());
  }
  LOG(1, "RALS") << "Using " << ansatz << " ansatz for decomposition"
                 << std::endl;
  computeOutgoingPi();

  out.set<Tensor<sisi4s::complex> *>("FactorOrbitals", PirR);
  out.set<Tensor<sisi4s::complex> *>("CoulombFactors", LambdaGR);
  if (isArgumentGiven("OutgoingFactorOrbitals")) {
    out.set<Tensor<sisi4s::complex> *>("OutgoingFactorOrbitals", PiqR);
  }
  composedGammaGqr = new Tensor<sisi4s::complex>(3,
                                                 GammaGqr->lens,
                                                 GammaGqr->sym,
                                                 *GammaGqr->wrld,
                                                 "composedGammaGqr",
                                                 GammaGqr->profile);
  if (isArgumentGiven("ComposedCoulombVertex")) {
    out.set<Tensor<sisi4s::complex> *>("ComposedCoulombVertex",
                                       composedGammaGqr);
  }

  double swampingThreshold(
      in.get<double>("swampingThreshold", DEFAULT_SWAMPING_THRESHOLD));
  double regularizationFriction(
      in.get<double>("regularizationFriction",
                     DEFAULT_REGULARIZATION_FRICTION));
  regularizationEstimator =
      new AlternatingLeastSquaresRegularizationEstimator(swampingThreshold,
                                                         regularizationFriction,
                                                         1);
  int64_t iterationsCount(0);
  int64_t maxIterationsCount(
      in.get<int64_t>("maxIterations", DEFAULT_MAX_ITERATIONS));
  double delta(in.get<double>("delta", DEFAULT_DELTA));
  Delta = std::numeric_limits<double>::infinity();
  while (iterationsCount < maxIterationsCount && Delta > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }

  // free up memroy
  if (!isArgumentGiven("OutgoingFactorOrbitals") && PiqR) { delete PiqR; }
  if (!isArgumentGiven("ComposedCoulombVertex") && composedGammaGqr) {
    delete composedGammaGqr;
  }
  if (regularizationEstimator) delete regularizationEstimator;
}

void CoulombVertexDecomposition::dryRun() {
  // NOTE that in the dry run GammaGai,... are local variables
  DryTensor<sisi4s::complex> *GammaGqr(
      in.get<DryTensor<sisi4s::complex> *>("CoulombVertex"));
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);

  // calculate decomposition rank
  rank = in.get<int64_t>("rankSize", DEFAULT_RANK_SIZE);
  // if rank is not given use rank factors (if they are not given use
  // rankFactors=2.0)
  if (rank == -1) {
    double rankFactor(in.get<double>("rankFactor", DEFAULT_RANK_FACTOR));
    rank = NG * rankFactor;
  }

  realFactorOrbitals =
      in.get<int64_t>("realFactorOrbitals", DEFAULT_REAL_FACTOR_ORBITALS);
  normalizedFactorOrbitals =
      in.get<int64_t>("normalizedFactorOrbitals",
                      DEFAULT_NORMALIZED_FACTOR_ORBITALS);
  LOG(0, "RALS") << "Tensor rank decomposition with rank NR=" << rank
                 << ", realFactorOrbitals=" << realFactorOrbitals
                 << ", normalizedFactorOrbitals=" << normalizedFactorOrbitals
                 << std::endl;
  LOG(1, "RALS") << "Decomposing Coulomb vertex with NG=" << NG << " Np=" << Np
                 << std::endl;

  if (isArgumentGiven("StartingFactorOrbitals")) {
    LOG(1, "RALS") << "Initial PirR=StartingPirR" << std::endl;
  } else {
    LOG(1, "RALS") << "Initial PirR=RandomTensor" << std::endl;
  }

  if (isArgumentGiven("StartingCoulombFactors")) {
    LOG(1, "RALS") << "Initial LambdaGR=StartingLambdaGR" << std::endl;
  } else {
    LOG(1, "RALS") << "Initial LambdaGR=RandomTensor" << std::endl;
  }

  // allocate factor tensors
  DryTensor<sisi4s::complex> *PiqR = new DryMatrix<complex>(Np, int(rank), NS);
  DryTensor<sisi4s::complex> *PirR = new DryMatrix<complex>(Np, int(rank), NS);
  DryTensor<sisi4s::complex> *LambdaGR =
      new DryMatrix<complex>(NG, int(rank), NS);
  out.set<DryTensor<sisi4s::complex> *>("FactorOrbitals", PirR);
  out.set<DryTensor<sisi4s::complex> *>("CoulombFactors", LambdaGR);

  std::string ansatz(in.get<std::string>("ansatz", HERMITIAN));
  if (ansatz == HERMITIAN) {
    LOG(1, "RALS") << "Using " << HERMITIAN << " ansatz for decomposition"
                   << std::endl;
  } else if (ansatz == SYMMETRIC) {
    LOG(1, "RALS") << "Using " << SYMMETRIC << " ansatz for decomposition"
                   << std::endl;
  } else if (ansatz == PSEUDO_INVERSE) {
    LOG(1, "RALS") << "Using " << PSEUDO_INVERSE << " ansatz for decomposition"
                   << std::endl;
  } else {
    std::stringstream stringStream;
    stringStream << "Unknown decomposition ansatz \"" << ansatz << "\"";
    throw new EXCEPTION(stringStream.str());
  }
  if (isArgumentGiven("OutgoingFactorOrbitals")) {
    out.set<DryTensor<sisi4s::complex> *>("OutgoingFactorOrbitals", PiqR);
  }

  DryTensor<sisi4s::complex> *composedGammaGqr(
      new DryTensor<sisi4s::complex>(*GammaGqr));
  if (isArgumentGiven("ComposedCoulombVertex")) {
    out.set<DryTensor<sisi4s::complex> *>("ComposedCoulombVertex",
                                          composedGammaGqr);
  }
  dryFit(GammaGqr, PiqR, PirR, LambdaGR, composedGammaGqr);
}

void CoulombVertexDecomposition::fit(int64_t const iterationsCount) {

  int fitFactorOrbitals(in.get<int64_t>("fitFactorOrbitals", 1));

  if (fitFactorOrbitals) { iterateQuadraticFactor(iterationsCount); }

  int fitCoulombFactors(in.get<int64_t>("fitCoulombFactors", 1));

  if (fitCoulombFactors) {
    fitRegularizedAlternatingLeastSquaresFactor(*GammaGqr,
                                                "Gqr",
                                                *PirR,
                                                'r',
                                                *PiqR,
                                                'q',
                                                *LambdaGR,
                                                'G',
                                                regularizationEstimator);
  }

  Delta = getDelta();
  LOG(0, "RALS") << "iteration=" << (iterationsCount + 1) << " Delta=" << Delta
                 << std::endl;
}

void CoulombVertexDecomposition::dryFit(
    DryTensor<sisi4s::complex> *GammaGqr,
    DryTensor<sisi4s::complex> *PiqR,
    DryTensor<sisi4s::complex> *PirR,
    DryTensor<sisi4s::complex> *LambdaGR,
    DryTensor<sisi4s::complex> *composedGammaGqr) {
  dryFitRegularizedAlternatingLeastSquaresFactor(*GammaGqr,
                                                 "Gqr",
                                                 *PiqR,
                                                 'q',
                                                 *LambdaGR,
                                                 'G',
                                                 *PirR,
                                                 'r');
  dryFitRegularizedAlternatingLeastSquaresFactor(*GammaGqr,
                                                 "Gqr",
                                                 *LambdaGR,
                                                 'G',
                                                 *PirR,
                                                 'r',
                                                 *PiqR,
                                                 'q');
  dryFitRegularizedAlternatingLeastSquaresFactor(*GammaGqr,
                                                 "Gqr",
                                                 *PirR,
                                                 'r',
                                                 *PiqR,
                                                 'q',
                                                 *LambdaGR,
                                                 'G');
  dryComposeCanonicalPolyadicDecompositionTensors(*LambdaGR,
                                                  *PiqR,
                                                  *PirR,
                                                  *composedGammaGqr);
}

void CoulombVertexDecomposition::normalizePi(Tensor<sisi4s::complex> &Pi) {
  Bivar_Function<complex> fDot(&sisi4s::dot<complex>);
  CTF::Vector<complex> norm(Pi.lens[0], *Pi.wrld);
  // norm["q"] = Pi["qR"] * conj(Pi["qR"])
  norm.contract(1.0, Pi, "qR", Pi, "qR", 0.0, "q", fDot);
  Tensor<sisi4s::complex> quotient(Pi);
  Univar_Function<complex> fSqrt(&sisi4s::sqrt<complex>);
  // quotient["qR"] = sqrt(norm["q"])
  quotient.sum(1.0, norm, "q", 0.0, "qR", fSqrt);
  Bivar_Function<complex> fDivide(&sisi4s::divide<complex>);
  // Pi["qR"] = Pi["qR"] / quotient["qR"]
  Pi.contract(1.0, Pi, "qR", quotient, "qR", 0.0, "qR", fDivide);
}

void CoulombVertexDecomposition::realizePi(Tensor<sisi4s::complex> &Pi) {
  Univar_Function<complex> fConj(&sisi4s::conj<complex>);
  Tensor<sisi4s::complex> conjX(Pi);
  // conjX["qR"] = conj(Pi["qR"])
  conjX.sum(1.0, Pi, "qR", 0.0, "qR", fConj);
  Pi["qR"] += conjX["qR"];
  Pi["qR"] *= 0.5;
}

void CoulombVertexDecomposition::iterateQuadraticFactor(int i) {
  // create a mixer
  std::string mixerName(in.get<std::string>("mixer", "LinearMixer"));
  PTR(Mixer<complex>) mixer(MixerFactory<complex>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  // initial guess is current Pi^q_R, composedGamma-Gamma is the residuum
  PTR(FockVector<complex>) Pi(NEW(FockVector<complex>,
                                  std::vector<PTR(Tensor<sisi4s::complex>)>(
                                      {NEW(Tensor<sisi4s::complex>, *PiqR)}),
                                  std::vector<std::string>({"qR"})));
  double initialDelta(getDelta());

  // Babylonian algorithm to solve quadratic form
  int maxSubIterationsCount(in.get<int64_t>("maxSubIterations", 8));
  int minSubIterationsCount(in.get<int64_t>("minSubIterations", 2));
  int j(0);
  double Delta(2 * initialDelta);
  while (j < minSubIterationsCount
         || (initialDelta < Delta && j < maxSubIterationsCount)) {
    // first, compute Pi_r^R by least-squares, fit keeping all other fixed
    fitAlternatingLeastSquaresFactor(*GammaGqr,
                                     "Gqr",
                                     *PiqR,
                                     'q',
                                     *LambdaGR,
                                     'G',
                                     *PirR,
                                     'r');
    if (realFactorOrbitals) realizePi(*PirR);
    if (normalizedFactorOrbitals) normalizePi(*PirR);
    // then get new Pi^q_R by conjugate transposition or pseudo inversion
    computeOutgoingPi();

    // compute difference to old Pi^q_R as residuum and append to mixer
    auto PiChange(Pi);
    Pi = NEW(FockVector<complex>,
             std::vector<PTR(Tensor<sisi4s::complex>)>(
                 {NEW(Tensor<sisi4s::complex>, *PiqR)}),
             std::vector<std::string>({"qR"}));
    *PiChange -= *Pi;
    mixer->append(Pi, PiChange);
    // get the mixer's best estimate for Pi^q_R
    Pi = NEW(FockVector<complex>, *mixer->get());
    // and write it to PiqR
    (*PiqR)["qR"] = (*Pi->get(0))["qR"];

    // compute fit error
    Delta = getDelta();
    if (writeSubIterations) {
      LOG(1, "Babylonian") << "|Pi^(" << (i + 1) << "," << (j + 1) << ")"
                           << "Pi*^(" << (i + 1) << "," << (j + 1) << ")"
                           << "Lambda^(n) - Gamma|=" << Delta << std::endl;
    }
    ++j;
  }
}

void CoulombVertexDecomposition::computeOutgoingPi() {
  std::string ansatz(in.get<std::string>("ansatz", HERMITIAN));

  if (ansatz == HERMITIAN) {
    Univar_Function<complex> fConj(&sisi4s::conj<complex>);
    PiqR->sum(1.0, *PirR, "qR", 0.0, "qR", fConj);
  } else if (ansatz == SYMMETRIC) {
    (*PiqR)["qR"] = (*PirR)["qR"];
  } else if (ansatz == PSEUDO_INVERSE) {
    (*PiqR)["qR"] = IterativePseudoInverse<complex>(*PirR).get()["Rq"];
  } else {
    std::stringstream stringStream;
    stringStream << "Unknown decomposition ansatz \"" << ansatz << "\"";
    throw new EXCEPTION(stringStream.str());
  }
}

double CoulombVertexDecomposition::getDelta() {
  composeCanonicalPolyadicDecompositionTensors(*LambdaGR,
                                               *PiqR,
                                               *PirR,
                                               *composedGammaGqr);
  (*composedGammaGqr)["Gqr"] -= (*GammaGqr)["Gqr"];
  double Delta(frobeniusNorm(*composedGammaGqr));
  (*composedGammaGqr)["Gqr"] += (*GammaGqr)["Gqr"];
  return Delta;
}

const std::string CoulombVertexDecomposition::SYMMETRIC("symmetric");
const std::string CoulombVertexDecomposition::HERMITIAN("hermitian");
const std::string CoulombVertexDecomposition::PSEUDO_INVERSE("pseudoInverse");
