#include <algorithms/PQRSCoulombIntegralsToVertex.hpp>

#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <DryTensor.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <util/CTF.hpp>
#include <util/MpiCommunicator.hpp>
#include <math/PseudoInverseSvd.hpp>
#include <extern/Lapack.hpp>


using namespace sisi4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(PQRSCoulombIntegralsToVertex);

PQRSCoulombIntegralsToVertex::PQRSCoulombIntegralsToVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PQRSCoulombIntegralsToVertex::~PQRSCoulombIntegralsToVertex() {
}

void PQRSCoulombIntegralsToVertex::run() {

  if (!isArgumentGiven("PQRSCoulombIntegrals")){
    throw new EXCEPTION("Need pqrs coulomb integrals");
  }
  if (!isArgumentGiven("HoleEigenEnergies")){
    throw new EXCEPTION("Need hole energies");
  }
  if (!isArgumentGiven("ParticleEigenEnergies")){
    throw new EXCEPTION("Need particle energies");
  }

  auto pqrs(getTensorArgument<>("PQRSCoulombIntegrals"));
  auto epsi(getTensorArgument<>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ParticleEigenEnergies"));
  auto No = epsi->lens[0];
  auto Nv = epsa->lens[0];
  auto Np = No + Nv;

  double *A, *C;
//  if (!Sisi4s::world->rank) {
    A = new double[Np*Np*Np*Np];
//  }
  LOG(0, "PQRS->vertex") << "we work with " << Np << " states\n";

  // ok we have to read the whole tensor in rank 0 only
  auto prqs(*pqrs);
  prqs["pqrs"] = (*pqrs)["prqs"];
  prqs.read_all(A);
  LOG(0, "PQRS") << "read all done\n";

  int n(Np*Np);
  int lwork(3*n);
  double *w = new double[n];
  double *work = new double[3*n];
  int info;

  dsyev_("V", "U", &n, A, &n, w, work, &lwork, &info);
  LOG(0,"INFO") << "Diagonalization sucessful: " << info << std::endl;
  size_t g(0);
  double thresh(1e-12);
  for (size_t m(0); m < n; m++) if (w[m] > thresh) g++;
  g = n;

  LOG(0, "PQRS") << g << " elements larger than " << thresh << std::endl;
  // this line is in principle useless because we work with a truncated matrix
  // but currently this is not working properly
  for (size_t i(0); i < n; i++){
    for (size_t m(0); m < n; m++){
      if (w[m] > thresh) A[i + m*n] = sqrt(w[m]) * A[i + m*n];
      else               A[i + m*n] = 0.0        * A[i + m*n];
    }
  }

  // write vertex info to yaml file
  std::string yamlout;
  yamlout += "version: 100\ntype: Tensor\nscalarType: Complex64\n";
  yamlout += "indices:\n  momentum:\n    type: halfGrid\n";
  yamlout += "dimensions:\n  - length:     ";
  yamlout += std::to_string(g);
  yamlout += "\n    type: AuxiliaryField\n  - length:     ";
  yamlout += std::to_string(Np);
  yamlout += "\n    type: State\n  - length:     ";
  yamlout += std::to_string(Np);
  yamlout += "\n    type: State\nelements:\n  type: IeeeBinaryFile\n";
  yamlout += "binary: 1\nunit: 1\nmetaData:\n  halfGrid: 1";

  std::ofstream vertyaml;
  vertyaml.open("CoulombVertex.yaml");
  vertyaml << yamlout;
  vertyaml.close();

  // write complex vertex data to file
  std::vector< std::complex<double> > out(n*g);
  for (size_t i(0); i < n; i++){
    size_t k(0);
    for (size_t j(n-g); j < n; j++){
      out[k++ + i*g] = A[i + j*n];
//      out[i + j*n] = { A[i + j*n], 0.0 };
    }
  }

  auto vertex = std::fstream("CoulombVertex.elements", std::ios::out | std::ios::binary);
  vertex.write((char*)&out[0], n*g*sizeof(std::complex<double>));
  vertex.close();


  // write eigenenergies to file
  std::vector<double> ea(Nv), ei(No);
  epsi->read_all(ei.data());
  epsa->read_all(ea.data());

  double fermiEnergy;
  if ( ei[No-1] < 0. && ea[0] > 0.) fermiEnergy = 0.0;
  else fermiEnergy = ( ea[0] + ei[No-1] ) / 2.0;


  char buf[50];
  yamlout.clear();
  yamlout += "version: 100\ntype: Tensor\nscalarType: Real64\n";
  yamlout += "dimensions:\n  - length:   ";
  yamlout += std::to_string(No+Nv);
  yamlout += "\n    type: State\nelements:\n  type: TextFile\n";
  yamlout += "unit: 1\nmetaData:\n  fermiEnergy:  ";
  sprintf(buf, "%16.14lf\n", fermiEnergy);
  yamlout += buf;
  yamlout += "  energies:";
  for (size_t ii = 0; ii < No; ii++){
    sprintf(buf, "\n    -  %16.14lf", ei[ii]);
    yamlout += buf;
  }
  for (size_t ii = 0; ii < Nv; ii++){
    sprintf(buf, "\n    -  %16.14lf", ea[ii]);
    yamlout += buf;
  }
  yamlout += "\n";
  std::ofstream eigyaml;
  eigyaml.open("EigenEnergies.yaml");
  eigyaml << yamlout;
  eigyaml.close();

  std::ofstream eigdat;
  eigdat.open("EigenEnergies.elements");
  yamlout.clear();
  for (size_t ii = 0; ii < No; ii++){
    sprintf(buf, "%16.14lf\n", ei[ii]);
    yamlout += buf;
  }
  for (size_t ii = 0; ii < Nv; ii++){
    sprintf(buf, "%16.14lf\n", ea[ii]);
    yamlout += buf;
  }
  eigdat << yamlout;
  eigdat.close();

}

