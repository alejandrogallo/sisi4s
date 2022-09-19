/*  BEGINNER's GUIDE TO LIBINT
 *  --------------------------
 *
 * When we load a basis set, we get a Basis Set object, which consists
 * of shells.
 * Every shell has some basis functions inside (contracted gaussians),
 * and these contracted gaussians are build out of uncontracted gaussian
 * orbitals.
 *
 * The integrals are also in general computed shell-wise, for instance
 * if s1 is the shell 1 of aug-cc-pvdz of Neon
 *
 *  S   8   1.00
 *    17880.0000000              0.0007380
 *     2683.0000000              0.0056770
 *      611.5000000              0.0288830
 *      173.5000000              0.1085400
 *       56.6400000              0.2909070
 *       20.4200000              0.4483240
 *        7.8100000              0.2580260
 *        1.6530000              0.0150630
 *
 * and s2 is the shell 4 of the same basis set
 *
 *  P   3   1.00
 *       28.3900000              0.0460870
 *        6.2700000              0.2401810
 *        1.6950000              0.5087440
 *
 *
 */
#include <string>
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <util/Libint.hpp>
#include <algorithms/HartreeFockFromGaussian.hpp>
#include <algorithms/OneBodyFromGaussian.hpp>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <Eigen/Eigenvalues>
#include <util/Tensor.hpp>
#include <util/Emitter.hpp>
#include <algorithms/HartreeFockFromCoulombIntegrals.hpp>
#define LOGGER(_l) LOG(_l, "HartreeFockFromGaussian")
#define IF_GIVEN(_l, ...) if (isArgumentGiven(_l)) { __VA_ARGS__ }

using namespace sisi4s;
ALGORITHM_REGISTRAR_DEFINITION(HartreeFockFromGaussian);
HartreeFockFromGaussian::HartreeFockFromGaussian(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
HartreeFockFromGaussian::~HartreeFockFromGaussian() {}

Tensor<double>
eigenToCtfMatrix(const Eigen::MatrixXd &m) {
  const int rank_m = int(Sisi4s::world->rank == 0); // rank mask
  int syms[] = {NS, NS}, lens[] = {(int)m.rows(), (int)m.cols()};
  std::vector<int64_t> indices(rank_m * m.rows() * m.cols());
  Tensor<double> t(2, lens, syms, *Sisi4s::world);
  std::iota(indices.begin(), indices.end(), 0);
  t.write(indices.size(), indices.data(), &m(0));
  return t;
}

double
getNuclearRepulsionEnergy(std::vector<libint2::Atom>& structure)
{
  unsigned int i, j;
  double enuc(0.0), r2(0.0);
  LOGGER(1) << "Calculating nuclear repulsion energy" << std::endl;
  for (i = 0; i < structure.size(); ++i) {
    for (j = i + 1; j < structure.size(); ++j) {
      r2 = 0.0;
      r2 += pow(structure[i].x - structure[j].x, 2);
      r2 += pow(structure[i].y - structure[j].y, 2);
      r2 += pow(structure[i].z - structure[j].z, 2);
      enuc += structure[i].atomic_number * structure[j].atomic_number
              / sqrt(r2);
    }
  }
  return enuc;
}

struct ShellInfo {
  size_t size, begin, end, l;
  inline ShellInfo(const libint2::BasisSet &shells, const size_t i) {
    size = shells[i].size();
    begin = shells.shell2bf()[i];
    end = size + begin;
    l = shells[i].contr[0].l;
  }
  inline size_t operator[](const size_t g) const { return g % size; }
};

Eigen::MatrixXd
getOneBodyIntegrals_(
  const libint2::BasisSet& shells,
  const libint2::Operator obtype,
  const std::vector<libint2::Atom>& atoms) {

  // Get number of basis set functions
  const size_t Np = shells.nbf();
  Eigen::MatrixXd result(Np, Np);

  // construct the overlap integrals engine
  libint2::Engine engine(obtype, shells.max_nprim(), shells.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const libint2::Atom& atom : atoms) {
      q.push_back(
        {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}}
      );
    }
    engine.set_params(q);
  }

  const auto& resultBuffer = engine.results();

  // loop over unique shell pairs p>=q
  for(size_t p = 0; p < shells.size(); ++p) {
  for(size_t q = 0; q <= p; ++q) {

    const ShellInfo P(shells, p), Q(shells, q);

    engine.compute(shells[p], shells[q]);

    // copy the result buffer of size (P.size*Q.size) to the
    // bufferMatrix, needed to do a block assignment to the Eigen matrix
    Eigen::Map<const Eigen::MatrixXd>
      bufferMatrix(resultBuffer[0], P.size, Q.size);
    result.block(P.begin, Q.begin, P.size, Q.size) = bufferMatrix;

    if (p != q) {
      // transpose and copy to the non-explicitly calculated block
      result.block(Q.begin, P.begin, Q.size, P.size) = bufferMatrix.transpose();
    }

  }  // q
  }  // p

  return result;
}

Eigen::MatrixXd
getTwoBodyFock(const libint2::BasisSet& shells, const Eigen::MatrixXd& D) {

  const size_t Np = shells.nbf();
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(Np, Np);

  // turn on the engine
  libint2::Engine engine(
    libint2::Operator::coulomb,
    shells.max_nprim(),
    shells.max_l(),
    0);

  // integrals[0] points to the target shell set after every call to
  // engine.compute
  const auto& integrals = engine.results();

  for(size_t _P = 0; _P < shells.size(); ++_P) {
  for(size_t _Q = 0; _Q < shells.size(); ++_Q) {
  for(size_t _R = 0; _R < shells.size(); ++_R) {
  for(size_t _S = 0; _S < shells.size(); ++_S) {
    // Fock matrix loop indices: _P, _Q
    // dens matrix loop indices: _R, _S
    // TODO: use simmetry for _R and _S
    // TODO: angular momentum checking
    const ShellInfo P(shells, _P),
                    Q(shells, _Q),
                    R(shells, _R),
                    S(shells, _S);

    // Coulomb contribution to the Fock matrix is from {_P,_Q,_R,_S} ints
    engine.compute(shells[_P], shells[_Q], shells[_R], shells[_S]);
    const double* Vprqs = integrals[0];
    // if all integrals screened out, skip to next quartet
    if (Vprqs == nullptr) continue;

    // Hartree part
    for(size_t p(P.begin), Iprqs = 0; p < P.end; ++p) {
    for(size_t q(Q.begin); q < Q.end; ++q) {
    for(size_t r(R.begin); r < R.end; ++r) {
    for(size_t s(S.begin); s < S.end; ++s, ++Iprqs) {
      G(p, q) += D(r, s) * 2.0 * Vprqs[Iprqs];
    } // s
    } // r
    } // q
    } // p

    // exchange contribution to the Fock matrix is from {_P,_R,_Q,_S} ints
    engine.compute(shells[_P], shells[_R], shells[_Q], shells[_S]);
    const double* Vprsq = integrals[0];

    // Exchange part
    for(size_t p(P.begin), Iprsq = 0; p < P.end; ++p) {
    for(size_t r(R.begin); r < R.end; ++r) {
    for(size_t q(Q.begin); q < Q.end; ++q) {
    for(size_t s(S.begin); s < S.end; ++s, ++Iprsq) {
      G(p, q) -= D(r, s) * Vprsq[Iprsq];
    } // s
    } // q
    } // r
    } // p

  } // _S
  } // _R
  } // _Q
  } // _P

  return G;
}


void HartreeFockFromGaussian::run() {

  std::vector<std::string> allArguments =
    { "xyzStructureFile"
    , "basisSet"
    , "CoreHamiltonian"
    , "energyDifference"
    , "maxIterations"
    , "numberOfElectrons"
    , "OrbitalCoefficients"
    , "OverlapMatrix"
    , "initialOrbitalCoefficients"
    , "HartreeFockEnergy"
    , "HoleEigenEnergies"
    , "ParticleEigenEnergies"
    };
  checkArgumentsOrDie(allArguments);

  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""))
                  , basisSet(getTextArgument("basisSet", "sto-3g"))
                  ;
  double electronicConvergence(getRealArgument("energyDifference", 1e-4));
  int numberOfElectrons(getIntegerArgument("numberOfElectrons", -1));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 16))
             , i
             , nBasisFunctions
             , No, Nv, Np
             ;

  LOGGER(1) << "maxIterations: " << maxIterations  << std::endl;
  LOGGER(1) << "ediff: " << electronicConvergence << std::endl;

  // Initialize libint
  LOGGER(1) << "libint2: " << LIBINT_VERSION << std::endl;
  LOGGER(1) << "MAX_AM: " << LIBINT_MAX_AM << std::endl;
  libint2::initialize();

  LOGGER(1) << "structure: " << xyzStructureFile << std::endl;
  std::ifstream structureFileStream(xyzStructureFile.c_str());
  std::vector<libint2::Atom> atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();

  if (numberOfElectrons == -1) {
    numberOfElectrons = 0;
    for (auto &atom: atoms) {
      std::cout << atom.atomic_number << std::endl;
      std::cout << atom.x << " " << atom.y << " " << atom.z << std::endl;
      numberOfElectrons += atom.atomic_number;
    }
  }

  LOGGER(1) << "natoms: " << atoms.size() << std::endl;
  LOGGER(1) << "nelec: " << numberOfElectrons << std::endl;

  // initializing basis set and outputing relevant information
  LOGGER(1) << "basis: " << basisSet << std::endl;
  libint2::BasisSet shells(basisSet, atoms);
  nBasisFunctions = shells.nbf();
  // restricted hartree fock
  No = numberOfElectrons/2;
  Nv = nBasisFunctions - No;
  Np = Nv + No;
  LOGGER(1) << "Initializing basis set.." << std::endl;
  LOGGER(1) << "No: " << No << std::endl;
  LOGGER(1) << "Nv: " << Nv << std::endl;
  LOGGER(1) << "#shells: " << shells.size() << std::endl;
  LOGGER(1) << "#functions: " << nBasisFunctions << std::endl;
  LOGGER(1) << "max_l: " << shells.max_l() << std::endl;
  LOGGER(1) << "Max functions in shell = " << shells.max_nprim() << std::endl;

  LOGGER(1) << "shell: l\tAOS\tCGS" << std::endl;
  i = 0;
  for (auto &shell: shells) {
    i++;
    LOGGER(1) << "shell " << i << ": "
      << shell.size() << "\t"
      << shell.nprim() << "\t"
      << shell.ncontr() << "\t"
      << std::endl;
  }


  LOGGER(1) << "Calculating overlaps" << std::endl;
  Eigen::MatrixXd S = getOneBodyIntegrals_(
    shells,
    libint2::Operator::overlap, atoms);

  LOGGER(1) << "Calculating kinetic integrals" << std::endl;
  Eigen::MatrixXd T = getOneBodyIntegrals_(
    shells,
    libint2::Operator::kinetic,
    atoms);
  LOGGER(1) << "T(" << T.rows() << "," << T.cols() << ")" << std::endl;


  LOGGER(1) << "Compute nuclear repulsion integrals" << std::endl;
  Eigen::MatrixXd V = getOneBodyIntegrals_(
    shells,
    libint2::Operator::nuclear,
    atoms);
  LOGGER(1) << "V(" << V.rows() << "," << V.cols() << ")" << std::endl;

  LOGGER(1) << "Calculating the core hamiltonian" << std::endl;
  Eigen::MatrixXd H = T + V;

  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);

  LOGGER(1)
    << "mem:Fock "
    << sizeof(double) * nBasisFunctions * nBasisFunctions / std::pow(2, 30.0)
    << " GB"
    << std::endl;


  Eigen::MatrixXd D(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd F(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd D_last = D;

  D *= 0.0;

  LOGGER(1) << "Setting initial density matrix" << std::endl;

  IF_GIVEN("initialOrbitalCoefficients",
    LOGGER(1) << "with initialOrbitalCoefficients" << std::endl;
    const auto ic_ctf(getTensorArgument<double>("initialOrbitalCoefficients"));
    const auto ic(toEigenMatrix(*ic_ctf));
    const auto C_occ(ic.leftCols(No));
    D = C_occ * C_occ.transpose();
  ) else {
    LOGGER(1) << "with whatever" << std::endl;
    D *= 0.0;
    for (unsigned i=0 ; i < Np ; i++) {
    for (unsigned j=i ; j < i+1 ; j++) {
      D(i,j) = 1;
    }
    }
  }

  unsigned int iter(0);
  double rmsd(0);
  double energyDifference(0);
  double ehf(0);
  double ehfLast(0);
  Eigen::MatrixXd eps, C;

  EMIT() << YAML::Key << "iterations"
         << YAML::Value << YAML::BeginSeq;

  do {

    ++iter;

    ehfLast = ehf;
    D_last = D;

    F = H;
    LOGGER(2) << "calculating fock matrix" << std::endl;
    F += getTwoBodyFock(shells, D);

    // solve F C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
      gen_eig_solver(F, S);
    eps = gen_eig_solver.eigenvalues();

    // C now has all eigenvectors from the shell of course
    C = gen_eig_solver.eigenvectors();

    // Computer the new density D = C(occ) . C(occ)T
    auto C_occ(C.leftCols(No));
    D = C_occ * C_occ.transpose();

    // Compute the current hartree fock energy for the current density matrix
    ehf = 0.0;
    for (unsigned i=0 ; i < nBasisFunctions; i++) {
      for (unsigned j=0 ; j < nBasisFunctions ; j++) {
        ehf += D(i,j) * (H(i,j) + F(i,j));
      }
    }

    energyDifference = ehf - ehfLast;
    rmsd = (D - D_last).norm();

    if (iter == 1) {
      LOG(1, "HartreeFockIt") << std::setprecision(15) << std::setw(10) <<
        "Iter"    << "\t" <<
        "E"       << "\t" <<
        "DeltaE"  << "\t" <<
        "RMS(D)"
        << std::endl;
    }

    LOG(1, "HartreeFockIt") <<
      iter              << "\t" <<
      ehf               << "\t" <<
      energyDifference  << "\t" <<
      rmsd              << "\t" <<
    std::endl;

    EMIT() <<
      YAML::BeginMap <<
        YAML::Key << "iteration" <<
        YAML::Value << iter <<
        YAML::Key << "energy" <<
        YAML::Value <<
          YAML::BeginMap <<
            YAML::Key << "value" <<
            YAML::Value << ehf <<
          YAML::EndMap <<
      YAML::EndMap;

  } while (
      ((fabs(energyDifference) > electronicConvergence) ||
       (fabs(rmsd) > electronicConvergence))        &&
      (iter < maxIterations)
      );

  EMIT() << YAML::EndSeq;
  for (unsigned int e; e<eps.size(); e++) {
    LOGGER(1) << "band " << e + 1 << " = " << eps(e,0) << std::endl;
  }

  double enuc = getNuclearRepulsionEnergy(atoms);
  LOGGER(1) << "energy=" << ehf + enuc << std::endl;
  LOGGER(1) << "nuclear energy: " << enuc << std::endl;

  int syms[] = {NS, NS}
    , o[]    = {(int)No}
    , v[]    = {(int)Nv}
    ;
  std::vector<int64_t> indices;
  const int rank_m = int(Sisi4s::world->rank == 0); // rank mask

  IF_GIVEN("CoreHamiltonian",
    LOGGER(1) << "Exporting CoreHamiltonian" << std::endl;
    allocatedTensorArgument<double>(
      "CoreHamiltonian", new Tensor<double>(eigenToCtfMatrix(H)));
  )

  IF_GIVEN("OrbitalCoefficients",
    LOGGER(1) << "Exporting OrbitalCoefficients" << std::endl;
    allocatedTensorArgument<double>(
      "OrbitalCoefficients", new Tensor<double>(eigenToCtfMatrix(C)));
  )

  IF_GIVEN("HoleEigenEnergies",
    // epsilon for holes
    indices.resize(rank_m * No);
    std::iota(indices.begin(), indices.end(), 0);
    auto epsi(new Tensor<double>(1, o, syms, *Sisi4s::world, "epsi"));
    epsi->write(indices.size(), indices.data(), &eps(0));
    allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  )


  IF_GIVEN("ParticleEigenEnergies",
    // epsilon for particles
    indices.resize(rank_m * Nv);
    std::iota(indices.begin(), indices.end(), 0);
    auto epsa(new Tensor<double>(1, v, syms, *Sisi4s::world, "epsa"));
    epsa->write(indices.size(), indices.data(), &eps(0) + No);
    allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);
  )

  IF_GIVEN("HartreeFockEnergy",
    CTF::Scalar<double> hfEnergy;
    hfEnergy[""] = (ehf + enuc);
    allocatedTensorArgument<double>("HartreeFockEnergy",
                                    new CTF::Scalar<double>(hfEnergy));
  )

  IF_GIVEN("OverlapMatrix",
    LOGGER(1) << "Exporting OverlapMatrix in the atomic basis" << std::endl;
    allocatedTensorArgument<double>(
      "OverlapMatrix", new Tensor<double>(eigenToCtfMatrix(S)));
  )

  libint2::finalize();

  EMIT() <<
    //--
    YAML::Key << "energy-convergence" <<
    YAML::Value << electronicConvergence <<
    //--
    YAML::Key << "energy" <<
    YAML::Value <<
      YAML::BeginMap <<
        YAML::Key << "value" <<
        YAML::Value << ehf + enuc <<
        YAML::Key << "nuclear" <<
        YAML::Value << enuc <<
      YAML::EndMap;

}
