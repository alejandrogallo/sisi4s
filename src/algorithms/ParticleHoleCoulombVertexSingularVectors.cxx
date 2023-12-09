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


DEFSPEC(ParticleHoleCoulombVertexSingularVectors,
        SPEC_IN({"reduction",
                 SPEC_VALUE_DEF("TODO: DOC", double, DEFAULT_REDUCTION)},
                {"fieldVariables",
                 SPEC_VALUE_DEF("TODO: DOC", int64_t, DEFAULT_FIELD_VARIABLES)},
                {"FullParticleHoleCoulombVertex",
                 SPEC_VARIN("TODO: DOC", Tensor<complex> *)}),
        SPEC_OUT());

IMPLEMENT_ALGORITHM(ParticleHoleCoulombVertexSingularVectors) {
  // read the particle-hole Coulomb vertex GammaGai
  // Its singular value decomposition is U.Sigma.W*
  // where W* is a matrix with the compound orbital index (a,i)
  Tensor<complex> *GammaGai(
      in.get<Tensor<complex> *>("FullParticleHoleCoulombVertex"));

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
  int NF(in.get<int64_t>("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(in.get<double>("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  // write singular vectors back to CTF
  CTF::Matrix<complex> U(USSUT);
  scaU->write(U);
  // slice singular vectors U corresponding to NF largest singular values S
  int start[] = {0, NG - NF}, end[] = {NG, NG};
  out.set<Tensor<complex> *>("ParticleHoleCoulombVertexSingularVectors",
                             new Tensor<complex>(U.slice(start, end)));

  // TODO: also write out the singular values
  delete[] SS;
}

void ParticleHoleCoulombVertexSingularVectors::dryRun() {
  // Read the Coulomb vertex GammaGai
  DryTensor<complex> *GammaGai(
      in.get<DryTensor<complex> *>("FullParticleHoleCoulombVertex"));

  DryTensor<complex> conjGammaGai(*GammaGai, SOURCE_LOCATION);

  int NG(GammaGai->lens[0]);
  // get number of field variables
  int NF(in.get<int64_t>("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(in.get<double>("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  out.set<DryTensor<complex> *>(
      "ParticleHoleCoulombVertexSingularVectors",
      new DryMatrix<complex>(NG, NF, NS, SOURCE_LOCATION));
}
