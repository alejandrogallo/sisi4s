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
#include <vector>
#include <libint2.hpp>
#include <algorithms/HartreeFock.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <ctf.hpp>
#include <numeric>      // std::iota

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(HartreeFock);
HartreeFock::HartreeFock(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
HartreeFock::~HartreeFock() {}


double
getNuclearRepulsionEnergy(std::vector<libint2::Atom>& structure)
{
  unsigned int i, j;
  double enuc(0.0), r2(0.0);
  LOG(1, "HartreeFock") << "Calculating nuclear repulsion energy" << std::endl;
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
getOneBodyIntegrals(
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

  // resultBuffer[0] points to the target shell set after every call to
  // engine.compute()
  const auto& resultBuffer = engine.results();

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
    const double* buf_1234 = resultBuffer[0];
    // if all integrals screened out, skip to next quartet
    if (buf_1234 == nullptr) continue;

    // Hartree part
    for(size_t p(P.begin), f1234 = 0; p < P.end; ++p) {
    for(size_t q(Q.begin); q < Q.end; ++q) {
    for(size_t r(R.begin); r < R.end; ++r) {
    for(size_t s(S.begin); s < S.end; ++s, ++f1234) {
      G(p, q) += D(r, s) * 2.0 * buf_1234[f1234];
    } // s
    } // r
    } // q
    } // p

    // exchange contribution to the Fock matrix is from {_P,_R,_Q,_S} ints
    engine.compute(shells[_P], shells[_R], shells[_Q], shells[_S]);
    const double* buf_1324 = resultBuffer[0];

    // Exchange part
    for(size_t p(P.begin), f1324 = 0; p < P.end; ++p) {
    for(size_t r(R.begin); r < R.end; ++r) {
    for(size_t q(Q.begin); q < Q.end; ++q) {
    for(size_t s(S.begin); s < S.end; ++s, ++f1324) {
      G(p, q) -= D(r, s) * buf_1324[f1324];
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



void HartreeFock::run() {

  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""));
  const std::string basisSet(getTextArgument("basisSet", "sto-3g"));
  unsigned int maxIterations(getIntegerArgument("maxIterations", 16));
  double electronicConvergence(getRealArgument("energyDifference", 1e-4));
  int numberOfElectrons(getIntegerArgument("numberOfElectrons", -1));
  unsigned int i;
  unsigned int nBasisFunctions;
  unsigned int No, Nv, Np;
  double enuc;

  LOG(1, "HartreeFock") << "maxIterations: " << maxIterations  << std::endl;
  LOG(1, "HartreeFock") << "ediff: " << electronicConvergence << std::endl;

  // Initialize libint
  LOG(1, "HartreeFock") << "libint2: " << LIBINT_VERSION << std::endl;
  LOG(1, "HartreeFock") << "MAX_AM: " << LIBINT_MAX_AM << std::endl;
  libint2::initialize();

  LOG(1, "HartreeFock") << "structure: " << xyzStructureFile << std::endl;
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

  LOG(1, "HartreeFock") << "natoms: " << atoms.size() << std::endl;
  LOG(1, "HartreeFock") << "nelec: " << numberOfElectrons << std::endl;

  // initializing basis set and outputing relevant information
  LOG(1, "HartreeFock") << "basis: " << basisSet << std::endl;
  libint2::BasisSet shells(basisSet, atoms);
  nBasisFunctions = shells.nbf();
  // restricted hartree fock
  No = numberOfElectrons/2;
  Nv = nBasisFunctions - No;
  Np = Nv + No;
  LOG(1, "HartreeFock") << "Initializing basis set.." << std::endl;
  LOG(1, "HartreeFock") << "No: " << No << std::endl;
  LOG(1, "HartreeFock") << "Nv: " << Nv << std::endl;
  LOG(1, "HartreeFock") << "#shells: " << shells.size() << std::endl;
  LOG(1, "HartreeFock") << "#functions: " << nBasisFunctions << std::endl;
  LOG(1, "HartreeFock") << "max_l: " << shells.max_l() << std::endl;
  LOG(1, "HartreeFock")
    << "Max functions in shell = " << shells.max_nprim() << std::endl;

  LOG(1, "HartreeFock") << "shell: l\tAOS\tCGS" << std::endl;
  i = 0;
  for (auto &shell: shells) {
    i++;
    LOG(1, "HartreeFock") << "shell " << i << ": "
      << shell.size() << "\t"
      << shell.nprim() << "\t"
      << shell.ncontr() << "\t"
      << std::endl;
  }

  enuc = getNuclearRepulsionEnergy(atoms);
  LOG(1, "HartreeFock") << "nuclear energy: " << enuc << std::endl;

  LOG(1, "HartreeFock") << "Calculating overlaps" << std::endl;
  Eigen::MatrixXd S = getOneBodyIntegrals(
    shells,
    libint2::Operator::overlap, atoms);

  LOG(1, "HartreeFock") << "Calculating kinetic integrals" << std::endl;
  Eigen::MatrixXd T = getOneBodyIntegrals(
    shells,
    libint2::Operator::kinetic,
    atoms);
  LOG(1, "HartreeFock")
    << "T(" << T.rows() << "," << T.cols() << ")" << std::endl;


  LOG(1, "HartreeFock")
    << "Compute nuclear repulsion integrals" << std::endl;
  Eigen::MatrixXd V = getOneBodyIntegrals(
    shells,
    libint2::Operator::nuclear,
    atoms);
  LOG(1, "HartreeFock")
    << "V(" << V.rows() << "," << V.cols() << ")" << std::endl;

  LOG(1, "HartreeFock") << "Calculating the core hamiltonian" << std::endl;
  Eigen::MatrixXd H = T + V;

  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);

  LOG(1, "HartreeFock")
    << "mem:Fock "
    << sizeof(double) * nBasisFunctions * nBasisFunctions / std::pow(2, 30.0)
    << " GB"
    << std::endl;


  Eigen::MatrixXd D(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd F(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd D_last = D;

  D *= 0.0;

  LOG(1, "HartreeFock") << "Setting initial density matrix" << std::endl;
  for ( i=0 ; i < nBasisFunctions ; i++) {
  for (unsigned j=i ; j < i+1 ; j++) {
    D(i,j) = 1;
  }
  }
  LOG(1, "HartreeFock") << "\tdone" << std::endl;

  unsigned int iter(0);
  double rmsd(0);
  double energyDifference(0);
  double ehf(0);
  double ehfLast(0);
  Eigen::MatrixXd eps, C;

  do {

    ++iter;

    ehfLast = ehf;
    D_last = D;

    F = H;
    LOG(2, "HartreeFock") << "calculating fock matrix" << std::endl;
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


  } while (
      ((fabs(energyDifference) > electronicConvergence) ||
       (fabs(rmsd) > electronicConvergence))        &&
      (iter < maxIterations)
      );

  for (unsigned int e; e<eps.size(); e++) {
    LOG(1, "HartreeFock")
      << "band " << e + 1 << " = " << eps(e,0) << std::endl;
  }

  LOG(1, "HartreeFock") << "energy=" << ehf + enuc << std::endl;

  int syms[] = {NS, NS};
  // export stuff
  int pp[] = {(int)Np, (int)Np};
  auto ctfCoefficients(new CTF::Tensor<double>(2, pp, syms, *Cc4s::world, "C"));
  std::vector<int64_t> indices;
  indices.resize(Np*Np);
  std::iota(indices.begin(), indices.end(), 0);
  ctfCoefficients->write(indices.size(), indices.data(), &C(0));

  int o[] = {(int)No}, v[] = {(int)Nv};
  // epsilon for holes
  indices.resize(No);
  std::iota(indices.begin(), indices.end(), 0);
  auto epsi(new CTF::Tensor<double>(1, o, syms, *Cc4s::world, "epsi"));
  epsi->write(indices.size(), indices.data(), &eps(0));
  // epsilon for particles
  indices.resize(Nv);
  std::iota(indices.begin(), indices.end(), 0);
  auto epsa(new CTF::Tensor<double>(1, v, syms, *Cc4s::world, "epsa"));
  epsa->write(indices.size(), indices.data(), &eps(0) + No);

  indices.resize(Np*Np*Np*Np);
  std::iota(indices.begin(), indices.end(), 0);

  allocatedTensorArgument<double>("OrbitalCoefficients", ctfCoefficients);
  allocatedTensorArgument<double>("HoleEigenEnergies", epsi);
  allocatedTensorArgument<double>("ParticleEigenEnergies", epsa);

  libint2::finalize();

}
