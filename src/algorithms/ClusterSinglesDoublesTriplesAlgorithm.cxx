#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <tcc/DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
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
  Mixer<F> *mixer( MixerFactory<F>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  // number of iterations for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );

  F e(0);
  for (int i(0); i < maxIterationsCount; ++i) {
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
  }

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "computing energy from given amplitudes" << std::endl;
    e = getEnergy(amplitudes);
  }

  storeAmplitudes(amplitudes, {"Singles", "Doubles", "Triples"});
  return e;
}


template <typename F>
PTR(FockVector<F>) ClusterSinglesDoublesTriplesAlgorithm::createAmplitudes(
  std::initializer_list<std::string> amplitudeNames,
  std::initializer_list<std::initializer_list<int>> amplitudeLens,
  std::initializer_list<std::string> amplitudeIndices
) {
  std::vector<PTR(CTF::Tensor<F>)> amplitudeTensors;
  auto lensIterator( amplitudeLens.begin() );
  for (auto name: amplitudeNames) {
    std::stringstream initialDataName;
    initialDataName << "initial" << name << "Amplitudes";
    if (isArgumentGiven(initialDataName.str())) {
      // use given amplitudes as initial amplitudes
      amplitudeTensors.push_back(
        NEW(CTF::Tensor<F>, *getTensorArgument<F>( initialDataName.str() ))
      );
    } else {
      // otherwise, use zeros as initial amplitudes
      std::vector<int> lens(*lensIterator);
      std::vector<int> syms(lens.size(), NS);
      amplitudeTensors.push_back(
        NEW(CTF::Tensor<F>,
          lens.size(), lens.data(), syms.data(), *Cc4s::world, "T"
        )
      );
    }
    ++lensIterator;
  }
  return NEW(FockVector<F>,
    amplitudeTensors.begin(), amplitudeTensors.end(),
    amplitudeIndices.begin(), amplitudeIndices.end()
  );
}


template <typename F>
void ClusterSinglesDoublesTriplesAlgorithm::storeAmplitudes(
  const PTR(const FockVector<F>) &amplitudes,
  std::initializer_list<std::string> names
) {
  int component(0);
  for (auto name: names) {
    if (isArgumentGiven(getDataName(name, "Amplitudes"))) {
      allocatedTensorArgument<F>(
        getDataName(name, "Amplitudes"),
        new Tensor<F>(*amplitudes->get(component))
      );
    }
    ++component;
  }
}


template <typename F>
void ClusterSinglesDoublesTriplesAlgorithm::estimateAmplitudesFromResiduum(
  const PTR(FockVector<F>) &residuum,
  const PTR(const FockVector<F>) &amplitudes
) {
  F levelShift( getRealArgument("levelShift", DEFAULT_LEVEL_SHIFT) );
  // level shifted division for left hand side
  class LevelShiftedDivision {
  public:
    LevelShiftedDivision(F shift_): shift(shift_) { }
    void operator()(F d, F &r) {
      r = -r / (d + shift);
    }
  protected:
    F shift;
  } levelShiftedDivision(levelShift);

  // apply level shifting on right hand side
  *residuum -= levelShift * *amplitudes;

  for (unsigned int i(0); i < residuum->componentTensors.size(); ++i) {
    auto R( residuum->get(i) );
    const char *indices( residuum->getIndices(i).c_str() );
    Tensor<F> D(false, *R);
    D.set_name("D");
    calculateExcitationEnergies(D, residuum->getIndices(i));

    // divide by -Delta to get new estimate for T
    CTF::Transform<F, F>(std::function<void(F, F &)>(levelShiftedDivision)) (
      D[indices], (*R)[indices]
    );
  }
}

// instantiate
template
void ClusterSinglesDoublesTriplesAlgorithm::estimateAmplitudesFromResiduum(
  const PTR(FockVector<double>) &residuum,
  const PTR(const FockVector<double>) &amplitudes
);

template
void ClusterSinglesDoublesTriplesAlgorithm::estimateAmplitudesFromResiduum(
  const PTR(FockVector<complex>) &residuum,
  const PTR(const FockVector<complex>) &amplitudes
);


template <typename F>
void ClusterSinglesDoublesTriplesAlgorithm::calculateExcitationEnergies(
  CTF::Tensor<F> &D, const std::string &indices
) {
  auto epsi(getTensorArgument<>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // convert to type F (either complex or double)
  Tensor<F> Fepsi(1, &epsi->lens[0], epsi->sym, *epsi->wrld, "Fepsi");
  // NOTE: just copies if both arguments are real
  toComplexTensor(*epsi, Fepsi);
  Tensor<F> Fepsa(1, &epsa->lens[0], epsa->sym, *epsa->wrld, "Fepsa");
  toComplexTensor(*epsa, Fepsa);

  // create excitation energy tensor
  int excitationLevel(D.order/2);
  for (int p(0); p < excitationLevel; ++p) {
    char epsaIndex[2] = {indices[p], 0};
    D[indices.c_str()] += Fepsa[epsaIndex];
    char epsiIndex[2] = {indices[excitationLevel+p], 0};
    D[indices.c_str()] -= Fepsi[epsiIndex];
  }
}

// instantiate
template
void ClusterSinglesDoublesTriplesAlgorithm::calculateExcitationEnergies(
  CTF::Tensor<double> &D, const std::string &indices
);
template
void ClusterSinglesDoublesTriplesAlgorithm::calculateExcitationEnergies(
  CTF::Tensor<complex> &D, const std::string &indices
);


template <typename F>
void ClusterSinglesDoublesTriplesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<F> &R
) {
  // Build D
  DryTensor<F> D(R, SOURCE_LOCATION);
}

// instantiate
template
void ClusterSinglesDoublesTriplesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<double> &R
);

template
void ClusterSinglesDoublesTriplesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<complex> &R
);

std::string ClusterSinglesDoublesTriplesAlgorithm::getCapitalizedAbbreviation() {
  std::string capitalizedAbbreviation(getAbbreviation());
  std::transform(
    capitalizedAbbreviation.begin(), capitalizedAbbreviation.end(),
    capitalizedAbbreviation.begin(), ::toupper
  );
  return capitalizedAbbreviation;
}


std::string ClusterSinglesDoublesTriplesAlgorithm::getDataName(
  const std::string &type, const std::string &data
) {
  std::stringstream dataName;
  dataName << getAbbreviation() << type << data;
  return dataName.str();
}
