#include <algorithms/ParticleHoleCoulombVertexSingularVectors.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <memory>

using namespace sisi4s;
using std::make_shared;
using std::shared_ptr;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombVertexSingularVectors);

ParticleHoleCoulombVertexSingularVectors::
    ParticleHoleCoulombVertexSingularVectors(
        std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

ParticleHoleCoulombVertexSingularVectors::
    ~ParticleHoleCoulombVertexSingularVectors() {}

void ParticleHoleCoulombVertexSingularVectors::run() {
  // read the particle-hole Coulomb vertex GammaGai
  // Its singular value decomposition is U.Sigma.W*
  // where W* is a matrix with the compound orbital index (a,i)
  Tensor<complex> *GammaGai(
      getTensorArgument<complex>("FullParticleHoleCoulombVertex"));

  // construct the conjugate of the Coulomb vertex GammaGai
  Tensor<complex> conjGammaGai(*GammaGai);
  conjugate(conjGammaGai);

  // contract away the orbitals away, leaving U.Sigma^2.U*
  int NG(GammaGai->lens[0]);
  CTF::Matrix<complex> USSUT(NG, NG, *GammaGai->wrld, "USSUT");
  LOG(1, "ParticleHoleCoulombVertexSingularVectors")
      << "Contracting over orbitals of " << GammaGai->get_name()
      << " to get U.Sigma^2.U*, with NG=" << NG << std::endl;
  USSUT["GH"] = conjGammaGai["Gai"] * (*GammaGai)["Hai"];

  // use ScaLapack routines to diagonalise the USSUT matrix, i.e. find U
  BlacsWorld world(USSUT.wrld->rank, USSUT.wrld->np);
  auto scaUSSUT(make_shared<ScaLapackMatrix<complex>>(USSUT, &world));
  auto scaU(make_shared<ScaLapackMatrix<complex>>(*scaUSSUT));
  ScaLapackHermitianEigenSystemDc<complex> eigenSystem(scaUSSUT, scaU);
  double *SS(new double[NG]);
  eigenSystem.solve(SS);

  // get number of field variables
  int NF(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  // write singular vectors back to CTF
  CTF::Matrix<complex> U(USSUT);
  scaU->write(U);
  // slice singular vectors U corresponding to NF largest singular values S
  int start[] = {0, NG - NF}, end[] = {NG, NG};
  allocatedTensorArgument<complex>("ParticleHoleCoulombVertexSingularVectors",
                                   new Tensor<complex>(U.slice(start, end)));

  // TODO: also write out the singular values
  delete[] SS;
}

void ParticleHoleCoulombVertexSingularVectors::dryRun() {
  // Read the Coulomb vertex GammaGai
  DryTensor<complex> *GammaGai(getTensorArgument<complex, DryTensor<complex>>(
      "FullParticleHoleCoulombVertex"));

  DryTensor<complex> conjGammaGai(*GammaGai, SOURCE_LOCATION);

  int NG(GammaGai->lens[0]);
  // get number of field variables
  int NF(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  allocatedTensorArgument<complex, DryTensor<complex>>(
      "ParticleHoleCoulombVertexSingularVectors",
      new DryMatrix<complex>(NG, NF, NS, SOURCE_LOCATION));
}
