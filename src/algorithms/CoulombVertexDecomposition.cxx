
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

ALGORITHM_REGISTRAR_DEFINITION(CoulombVertexDecomposition);

CoulombVertexDecomposition::CoulombVertexDecomposition(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList)
    , composedGammaGqr(nullptr)
    , PiqR(nullptr)
    , regularizationEstimator(nullptr) {}

CoulombVertexDecomposition::~CoulombVertexDecomposition() {
  if (!isArgumentGiven("OutgoingFactorOrbitals") && PiqR) { delete PiqR; }
  if (!isArgumentGiven("ComposedCoulombVertex") && composedGammaGqr) {
    delete composedGammaGqr;
  }
  if (regularizationEstimator) delete regularizationEstimator;
}

void CoulombVertexDecomposition::run() {
  GammaGqr = getTensorArgument<complex>("CoulombVertex");
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);

  // calculate decomposition rank
  rank = getIntegerArgument("rankSize", DEFAULT_RANK_SIZE);
  // if rank is not given use rank factors (if they are not given use
  // rankFactors=3.0)
  if (rank == -1) {
    double rankFactor(getRealArgument("rankFactor", DEFAULT_RANK_FACTOR));
    rank = NG * rankFactor;
  }

  realFactorOrbitals =
      getIntegerArgument("realFactorOrbitals", DEFAULT_REAL_FACTOR_ORBITALS);
  normalizedFactorOrbitals =
      getIntegerArgument("normalizedFactorOrbitals",
                         DEFAULT_NORMALIZED_FACTOR_ORBITALS);
  LOG(0, "RALS") << "Tensor rank decomposition with rank NR=" << rank
                 << ", realFactorOrbitals=" << realFactorOrbitals
                 << ", normalizedFactorOrbitals=" << normalizedFactorOrbitals
                 << std::endl;
  LOG(1, "RALS") << "Decomposing Coulomb vertex " << GammaGqr->get_name()
                 << " with NG=" << NG << ", Np=" << Np << std::endl;

  writeSubIterations =
      getIntegerArgument("writeSubIterations", DEFAULT_WRITE_SUB_ITERATIONS);

  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);

  // allocate factor tensors
  if (isArgumentGiven("StartingFactorOrbitals")) {
    Tensor<complex> *PirRTensor(
        getTensorArgument<complex>("StartingFactorOrbitals"));
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
    Tensor<complex> *LambdaGRTensor(
        getTensorArgument<complex>("StartingCoulombFactors"));
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
  std::string ansatz(getTextArgument("ansatz", HERMITIAN));
  if (ansatz != HERMITIAN && ansatz != SYMMETRIC && ansatz != PSEUDO_INVERSE) {
    std::stringstream stringStream;
    stringStream << "Unknown decomposition ansatz \"" << ansatz << "\"";
    throw new EXCEPTION(stringStream.str());
  }
  LOG(1, "RALS") << "Using " << ansatz << " ansatz for decomposition"
                 << std::endl;
  computeOutgoingPi();

  allocatedTensorArgument<complex>("FactorOrbitals", PirR);
  allocatedTensorArgument<complex>("CoulombFactors", LambdaGR);
  if (isArgumentGiven("OutgoingFactorOrbitals")) {
    allocatedTensorArgument<complex>("OutgoingFactorOrbitals", PiqR);
  }
  composedGammaGqr = new Tensor<complex>(3,
                                         GammaGqr->lens,
                                         GammaGqr->sym,
                                         *GammaGqr->wrld,
                                         "composedGammaGqr",
                                         GammaGqr->profile);
  if (isArgumentGiven("ComposedCoulombVertex")) {
    allocatedTensorArgument<complex>("ComposedCoulombVertex", composedGammaGqr);
  }

  double swampingThreshold(
      getRealArgument("swampingThreshold", DEFAULT_SWAMPING_THRESHOLD));
  double regularizationFriction(
      getRealArgument("regularizationFriction",
                      DEFAULT_REGULARIZATION_FRICTION));
  regularizationEstimator =
      new AlternatingLeastSquaresRegularizationEstimator(swampingThreshold,
                                                         regularizationFriction,
                                                         1);
  int64_t iterationsCount(0);
  int64_t maxIterationsCount(
      getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  double delta(getRealArgument("delta", DEFAULT_DELTA));
  Delta = std::numeric_limits<double>::infinity();
  while (iterationsCount < maxIterationsCount && Delta > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }
}

void CoulombVertexDecomposition::dryRun() {
  // NOTE that in the dry run GammaGai,... are local variables
  DryTensor<complex> *GammaGqr(
      getTensorArgument<complex, DryTensor<complex>>("CoulombVertex"));
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);

  // calculate decomposition rank
  rank = getIntegerArgument("rankSize", DEFAULT_RANK_SIZE);
  // if rank is not given use rank factors (if they are not given use
  // rankFactors=2.0)
  if (rank == -1) {
    double rankFactor(getRealArgument("rankFactor", DEFAULT_RANK_FACTOR));
    rank = NG * rankFactor;
  }

  realFactorOrbitals =
      getIntegerArgument("realFactorOrbitals", DEFAULT_REAL_FACTOR_ORBITALS);
  normalizedFactorOrbitals =
      getIntegerArgument("normalizedFactorOrbitals",
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
  DryTensor<complex> *PiqR = new DryMatrix<complex>(Np, int(rank), NS);
  DryTensor<complex> *PirR = new DryMatrix<complex>(Np, int(rank), NS);
  DryTensor<complex> *LambdaGR = new DryMatrix<complex>(NG, int(rank), NS);
  allocatedTensorArgument<complex, DryTensor<complex>>("FactorOrbitals", PirR);
  allocatedTensorArgument<complex, DryTensor<complex>>("CoulombFactors",
                                                       LambdaGR);

  std::string ansatz(getTextArgument("ansatz", HERMITIAN));
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
    allocatedTensorArgument<complex, DryTensor<complex>>(
        "OutgoingFactorOrbitals",
        PiqR);
  }

  DryTensor<complex> *composedGammaGqr(new DryTensor<complex>(*GammaGqr));
  if (isArgumentGiven("ComposedCoulombVertex")) {
    allocatedTensorArgument<complex, DryTensor<complex>>(
        "ComposedCoulombVertex",
        composedGammaGqr);
  }
  dryFit(GammaGqr, PiqR, PirR, LambdaGR, composedGammaGqr);
}

void CoulombVertexDecomposition::fit(int64_t const iterationsCount) {

  int fitFactorOrbitals(getIntegerArgument("fitFactorOrbitals", 1));

  if (fitFactorOrbitals) { iterateQuadraticFactor(iterationsCount); }

  int fitCoulombFactors(getIntegerArgument("fitCoulombFactors", 1));

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

void CoulombVertexDecomposition::dryFit(DryTensor<complex> *GammaGqr,
                                        DryTensor<complex> *PiqR,
                                        DryTensor<complex> *PirR,
                                        DryTensor<complex> *LambdaGR,
                                        DryTensor<complex> *composedGammaGqr) {
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

void CoulombVertexDecomposition::normalizePi(Tensor<complex> &Pi) {
  Bivar_Function<complex> fDot(&sisi4s::dot<complex>);
  CTF::Vector<complex> norm(Pi.lens[0], *Pi.wrld);
  // norm["q"] = Pi["qR"] * conj(Pi["qR"])
  norm.contract(1.0, Pi, "qR", Pi, "qR", 0.0, "q", fDot);
  Tensor<complex> quotient(Pi);
  Univar_Function<complex> fSqrt(&sisi4s::sqrt<complex>);
  // quotient["qR"] = sqrt(norm["q"])
  quotient.sum(1.0, norm, "q", 0.0, "qR", fSqrt);
  Bivar_Function<complex> fDivide(&sisi4s::divide<complex>);
  // Pi["qR"] = Pi["qR"] / quotient["qR"]
  Pi.contract(1.0, Pi, "qR", quotient, "qR", 0.0, "qR", fDivide);
}

void CoulombVertexDecomposition::realizePi(Tensor<complex> &Pi) {
  Univar_Function<complex> fConj(&sisi4s::conj<complex>);
  Tensor<complex> conjX(Pi);
  // conjX["qR"] = conj(Pi["qR"])
  conjX.sum(1.0, Pi, "qR", 0.0, "qR", fConj);
  Pi["qR"] += conjX["qR"];
  Pi["qR"] *= 0.5;
}

void CoulombVertexDecomposition::iterateQuadraticFactor(int i) {
  // create a mixer
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  PTR(Mixer<complex>) mixer(MixerFactory<complex>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  // initial guess is current Pi^q_R, composedGamma-Gamma is the residuum
  PTR(FockVector<complex>) Pi(
      NEW(FockVector<complex>,
          std::vector<PTR(Tensor<complex>)>({NEW(Tensor<complex>, *PiqR)}),
          std::vector<std::string>({"qR"})));
  double initialDelta(getDelta());

  // Babylonian algorithm to solve quadratic form
  int maxSubIterationsCount(getIntegerArgument("maxSubIterations", 8));
  int minSubIterationsCount(getIntegerArgument("minSubIterations", 2));
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
             std::vector<PTR(Tensor<complex>)>({NEW(Tensor<complex>, *PiqR)}),
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
  std::string ansatz(getTextArgument("ansatz", HERMITIAN));

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
