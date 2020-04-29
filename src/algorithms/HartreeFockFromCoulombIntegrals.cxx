#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <string>
#include <vector>
#include <algorithms/HartreeFockFromCoulombIntegrals.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#define IF_GIVEN(_l, ...) if (isArgumentGiven(_l)) { __VA_ARGS__ }
#define LOGGER(_l) LOG(_l, "HartreeFockFromCoulombIntegrals")
#define LOGGER_IT(_l) LOG(_l, "HartreeFockFromCoulombIntegralsIt")


using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(HartreeFockFromCoulombIntegrals);

Eigen::MatrixXd
cc4s::toEigenMatrix(CTF::Tensor<double> &ctf) {
  Eigen::MatrixXd result(ctf.lens[0], ctf.lens[1]);

  const size_t size(ctf.lens[0] * ctf.lens[1]);
  std::vector<int64_t> indices(size);
  std::vector<double> values(size);

  // make indices a range from 0 to indices.size()-1
  std::iota(indices.begin(), indices.end(), 0);

  ctf.read(indices.size(), indices.data(), values.data());

  for (size_t i(0); i<values.size(); i++) {
    result(i) = values[i];
  }

  return result;
}

CTF::Tensor<double>
toCtfMatrix(const Eigen::MatrixXd &m) {
  int syms[] = {NS, NS}, lens[] = {(int)m.rows(), (int)m.cols()};
  LOGGER(2) << "converting into ctf vector" << std::endl;
  const int64_t size(m.rows() * m.cols())
              , r(Cc4s::world->rank)
              , np(Cc4s::world->np)
              , chunks(size / np)
              , start(r * chunks)
              , end((np-1) == r ? size : (r+1)*chunks);
              ;
  std::vector<int64_t> indices(end - start);
  CTF::Tensor<double> t(2, lens, syms, *Cc4s::world);

  //  make indices a range from 0 to indices.size()-1
  std::iota(indices.begin(), indices.end(), start);
  t.write(indices.size(), indices.data(), &m(0) + start);

  return t;
}

Eigen::MatrixXd
getFockMatrix(const Eigen::MatrixXd &d, CTF::Tensor<double> &V) {
  auto D(toCtfMatrix(d));
  CTF::Tensor<double> F(2, D.lens, D.sym, *Cc4s::world);
  LOGGER(1) << "Hartree" << std::endl;
  F["pq"]  = ( 2.0) * D["kl"] * V["kplq"];
  LOGGER(1) << "Fock" << std::endl;
  F["pq"] += (-1.0) * D["kl"] * V["kpql"];
  return cc4s::toEigenMatrix(F);
}

void HartreeFockFromCoulombIntegrals::run() {
  checkArgumentsOrDie( { "h"
                       , "CoulombIntegrals"
                       , "No"
                       , "Nv"
                       , "maxIterations"
                       , "initialOrbitalCoefficients"
                       , "energyDifference"
                       , "OverlapMatrix"
                       , "OrbitalCoefficients"
                       , "HartreeFockEnergy"
                       , "HoleEigenEnergies"
                       , "ParticleEigenEnergies"
                       } );

  const auto ctfH(getTensorArgument<double>("h"));
  const auto H(cc4s::toEigenMatrix(*ctfH));
  const auto V(getTensorArgument<double>("CoulombIntegrals"));
  const unsigned int No(getIntegerArgument("No"));
  const size_t Nv(getIntegerArgument("Nv", V->lens[0] - No));
  const size_t Np(No+Nv);
  const unsigned int maxIterations(getIntegerArgument("maxIterations", 16));
  const double electronicConvergence(getRealArgument("energyDifference", 1e-4));

  LOGGER(1) << "maxIterations: " << maxIterations  << std::endl;
  LOGGER(1) << "ediff: " << electronicConvergence << std::endl;
  LOGGER(1) << "No: " << No << std::endl;
  LOGGER(1) << "Nv: " << Nv << std::endl;
  LOGGER(1) << "Calculating overlaps" << std::endl;

  Eigen::MatrixXd S(Eigen::MatrixXd::Identity(Np, Np));
  if (isArgumentGiven("OverlapMatrix")) {
    LOGGER(1) << "Setting provided OverlapMatrix" << std::endl;
    S = cc4s::toEigenMatrix(*getTensorArgument<double>("OverlapMatrix"));
  }

  LOGGER(1) << "mem:Fock "
            << sizeof(double) * Np * Np / std::pow(2, 30.0)
            << " GB"
            << std::endl;

  MatrixColumnMajor D(Np, Np), F(Np, Np), fockMatrix(Np, Np), D_last(D);


  LOGGER(1) << "Setting initial density matrix" << std::endl;

  IF_GIVEN("initialOrbitalCoefficients",
    LOGGER(1) << "with initialOrbitalCoefficients" << std::endl;
    const auto ic_ctf(getTensorArgument<double>("initialOrbitalCoefficients"));
    const auto ic(cc4s::toEigenMatrix(*ic_ctf));
    const auto C_occ(ic.leftCols(No));
    D = C_occ * C_occ.transpose();
  ) else {
    LOGGER(1) << "diagonalizing overlap matrix" << std::endl;
    Eigen::SelfAdjointEigenSolver<MatrixColumnMajor> solver(S);
    auto C(solver.eigenvectors());
    auto C_occ(C.leftCols(No));
    LOGGER(1) << "with eigenvectors" << std::endl;
    D = C_occ * C_occ.transpose();
  }

  unsigned int iter(0);
  double rmsd(0);
  double energyDifference(0);
  double ehf(0);
  double ehfLast(0);
  MatrixColumnMajor eps, C;

  const auto updateHamiltonian = [&] { F = H;
                                       fockMatrix = getFockMatrix(D, *V);
                                       F += fockMatrix;
                                     };

  const auto updateEnergy = [&] { ehf = 0.0;
                                  for (size_t i=0 ; i < Np; i++)
                                  for (size_t j=0 ; j < Np ; j++)
                                    ehf += D(i,j) * (H(i,j) + F(i,j));
                                  // update energy difference
                                  energyDifference = ehf - ehfLast;
                                };

  const auto updateDensity = [&] {
    LOGGER(1) << "Diagonalize" << std::endl;
    // solve F C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixColumnMajor>
      gen_eig_solver(F, S);
    eps = gen_eig_solver.eigenvalues();

    // C now has all eigenvectors from the shell of course
    C = gen_eig_solver.eigenvectors();

    // Computer the new density D = C(occ) . C(occ)T
    auto C_occ(C.leftCols(No));
    D = C_occ * C_occ.transpose();

    rmsd = (D - D_last).norm();
  };

  LOGGER(2) << "calculating initial fock matrix" << std::endl;
  updateHamiltonian();
  updateEnergy();
  LOGGER(1) << "initial guess energy = " << ehf << std::endl;

  do {

    ++iter;

    ehfLast = ehf;
    D_last = D;

    updateHamiltonian();
    updateDensity();
    updateEnergy();


    if (iter == 1) {
      LOGGER(1) << std::setprecision(16) << std::setw(10) <<
        "Iter"    << "\t" <<
        "E"       << "\t" <<
        "DeltaE"  << "\t" <<
        "RMS(D)"
        << std::endl;
    }

    LOGGER_IT(1) <<
      iter              << "\t" <<
      ehf               << "\t" <<
      energyDifference  << "\t" <<
      rmsd              << "\t" <<
    std::endl;


  } while ( (  (fabs(energyDifference) > electronicConvergence)
            || (fabs(rmsd) > electronicConvergence)
            )
          && (iter < maxIterations)
          );

  for (unsigned int e; e<eps.size(); e++) {
    LOGGER(1) << "band " << e + 1 << " = " << eps(e,0) << std::endl;
  }

  LOGGER(1) << "energy=" << ehf << std::endl;

  // export stuff
  const int rank_m = int(Cc4s::world->rank == 0); // rank mask
  int syms[] = {NS, NS};
  int o[] = {(int)No}, v[] = {(int)Nv};
  std::vector<int64_t> indices;

  IF_GIVEN("FockMatrix",
    int pp[] = {(int)Np, (int)Np};
    auto f(new CTF::Tensor<double>(2, pp, syms, *Cc4s::world, "F"));
    indices.resize(rank_m * Np * Np);
    std::iota(indices.begin(), indices.end(), 0);
    f->write(indices.size(), indices.data(), &fockMatrix(0));
    allocatedTensorArgument<double>("FockMatrix", f);
  )

  IF_GIVEN("OrbitalCoefficients",
    int pp[] = {(int)Np, (int)Np};
    auto ctfcs(new CTF::Tensor<double>(2, pp, syms, *Cc4s::world, "C"));
    indices.resize(rank_m * Np * Np);
    std::iota(indices.begin(), indices.end(), 0);
    ctfcs->write(indices.size(), indices.data(), &C(0));
    allocatedTensorArgument<double>("OrbitalCoefficients", ctfcs);
  )

  IF_GIVEN("HoleEigenEnergies",
    // epsilon for holes
    indices.resize(rank_m * No);
    std::iota(indices.begin(), indices.end(), 0);
    auto epsi(new CTF::Tensor<double>(1, o, syms, *Cc4s::world, "epsi"));
    epsi->write(indices.size(), indices.data(), &eps(0));
    allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  )

  IF_GIVEN("ParticleEigenEnergies",
    // epsilon for particles
    indices.resize(rank_m * Nv);
    std::iota(indices.begin(), indices.end(), 0);
    auto epsa(new CTF::Tensor<double>(1, v, syms, *Cc4s::world, "epsa"));
    epsa->write(indices.size(), indices.data(), &eps(0) + No);
    allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);
  )

  IF_GIVEN("HartreeFockEnergy",
    auto hfEnergy = new CTF::Scalar<double>();
    (*hfEnergy)[""] = (ehf);
    allocatedTensorArgument<double>("HartreeFockEnergy", hfEnergy);
  )

}
