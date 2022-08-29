#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <tcc/DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/CTF.hpp>
#include <Options.hpp>
#include <Cc4s.hpp>

#include <initializer_list>

using namespace CTF;
using namespace cc4s;


ClusterSinglesDoublesTriplesAlgorithm::~ClusterSinglesDoublesTriplesAlgorithm() {
}

ClusterSinglesDoublesTriplesAlgorithm::ClusterSinglesDoublesTriplesAlgorithm(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}


void ClusterSinglesDoublesTriplesAlgorithm::run() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = run<double>();
  } else {
    e = std::real( run<complex>() );
  }
  setRealArgument(getDataName("", "Energy"), e);
}

template <typename F>
F ClusterSinglesDoublesTriplesAlgorithm::run() {
  int Nv(getTensorArgument<>("ParticleEigenEnergies")->lens[0]);
  int No(getTensorArgument<>("HoleEigenEnergies")->lens[0]);

  PTR(const FockVector<F>) amplitudes(
    createAmplitudes<F>(
      {"Singles", "Doubles", "Triples"},
      {{Nv,No}, {Nv,Nv,No,No}, {Nv,Nv,Nv,No,No,No}},
      {"ai", "abij", "abcijk"}
    )
  );

  // create a mixer, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  PTR(Mixer<F>) mixer( MixerFactory<F>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  // number of iterations for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );

  F amplitudesConvergence(
    getRealArgument("amplitudesConvergence", DEFAULT_AMPLITUDES_CONVERGENCE)
  );
  F energyConvergence(
    getRealArgument("energyConvergence", DEFAULT_ENERGY_CONVERGENCE)
  );

  F e(0), previousE(0);
  int i(0);
  for (; i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    // call the getResiduum of the actual algorithm,
    // which will be specified by inheriting classes
    auto estimatedAmplitudes( getResiduum(i, amplitudes) );
    estimateAmplitudesFromResiduum(estimatedAmplitudes, amplitudes);
    auto amplitudesChange( NEW(FockVector<F>, *estimatedAmplitudes) );
    *amplitudesChange -= *amplitudes;
    mixer->append(estimatedAmplitudes, amplitudesChange);
    // get mixer's best guess for amplitudes
    amplitudes = mixer->get();
    e = getEnergy(amplitudes);
    if (
      std::abs((e-previousE)/e) < std::abs(energyConvergence) &&
      std::abs(
        amplitudesChange->dot(*amplitudesChange) / amplitudes->dot(*amplitudes)
      ) < std::abs(amplitudesConvergence * amplitudesConvergence)
    ) break;
    previousE = e;
  }

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "computing energy from given amplitudes" << std::endl;
    e = getEnergy(amplitudes);
  } else if (i == maxIterationsCount) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "WARNING: energy or amplitudes convergence not reached." << std::endl;
  }

  storeAmplitudes(amplitudes, {"Singles", "Doubles", "Triples"});
  return e;
}
