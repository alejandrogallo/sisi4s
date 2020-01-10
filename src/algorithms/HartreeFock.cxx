#include <string>
#include <vector>
#include <libint2.hpp>
#include <algorithms/HartreeFock.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Eigen/Eigenvalues>

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


Eigen::MatrixXd
getOneBodyIntegrals(
  const libint2::BasisSet shells,
  libint2::Operator obtype,
  const std::vector<libint2::Atom>& atoms
  )
{

  // Get number of basis set functions
  int n = shells.nbf();
  Eigen::MatrixXd result(n,n);

  // construct the overlap integrals engine
  libint2::Engine engine(obtype, shells.max_nprim(), shells.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const libint2::Atom& atom : atoms) {
      q.push_back(
          {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}}
          );
    }
    engine.set_params(q);
  }

  const std::vector<size_t> shell2bf = shells.shell2bf();

  // buf[0] points to the target shell set after every call  to
  // engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for(size_t s1=0; s1!=shells.size(); ++s1) {

    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(size_t s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
      //block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;
}


Eigen::MatrixXd
getTwoBodyFock(
  const libint2::BasisSet shells,
  const Eigen::MatrixXd& D
  ) {

  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = shells.nbf();
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n,n);

  // construct the electron repulsion integrals engine
  libint2::Engine engine(
    libint2::Operator::coulomb,
    shells.max_nprim(),
    shells.max_l(),
    0
  );

  const std::vector<size_t> shell2bf = shells.shell2bf();

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
  for(size_t s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(size_t s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(size_t s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(size_t s4=0; s4!=shells.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr)
            continue; // if all integrals screened out, skip to next quartet

          // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(size_t f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(size_t f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(size_t f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(size_t f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * buf_1234[f1234];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
          const auto* buf_1324 = buf[0];

          for(size_t f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(size_t f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(size_t f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(size_t f4=0; f4!=n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) -= D(bf3,bf4) * buf_1324[f1324];
                }
              }
            }
          }

        }
      }
    }
  }

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
  unsigned int No;
  double enuc;

  LOG(1, "HartreeFock") << "maxIterations: " << maxIterations  << std::endl;
  LOG(1, "HartreeFock") << "ediff: " << electronicConvergence << std::endl;
  LOG(1, "HartreeFock") << "Basis set: " << basisSet << std::endl;

  // Initialize libint
  LOG(1, "HartreeFock") << "libint2: " << LIBINT_VERSION << std::endl;
  LOG(1, "HartreeFock") << "MAX_AM: " << LIBINT_MAX_AM << std::endl;
  libint2::initialize();

  LOG(1, "HartreeFock") << "structure: " << xyzStructureFile << std::endl;
  std::ifstream structureFileStream(xyzStructureFile.c_str());
  std::vector<libint2::Atom> atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();

  libint2::BasisSet shells(basisSet, atoms);
  nBasisFunctions = shells.nbf();

  if (numberOfElectrons == -1) {
    numberOfElectrons = 0;
    for (auto &atom: atoms) {
      std::cout << atom.atomic_number << std::endl;
      std::cout << atom.x << " " << atom.y << " " << atom.z << std::endl;
      numberOfElectrons += atom.atomic_number;
    }
  }

  // restricted hartree fock
  No = numberOfElectrons/2;
  LOG(1, "HartreeFock") << "natoms: " << atoms.size() << std::endl;
  LOG(1, "HartreeFock") << "nelec: " << numberOfElectrons << std::endl;
  LOG(1, "HartreeFock") << "No: " << No << std::endl;

  enuc = getNuclearRepulsionEnergy(atoms);
  LOG(1, "HartreeFock") << "nuclear energy: " << enuc << std::endl;

  LOG(1, "HartreeFock") << "Initializing basis set.." << std::endl;

  LOG(1, "HartreeFock") << "#basis: " << nBasisFunctions << std::endl;
  LOG(1, "HartreeFock")
    << "Maximum primitives in shells = " << shells.max_nprim()
    << std::endl;
  LOG(1, "HartreeFock") << "max_l: " << shells.max_l() << std::endl;

  LOG(1, "HartreeFock") << "Calculating overlaps" << std::endl;
  Eigen::MatrixXd S = getOneBodyIntegrals(
    shells,
    libint2::Operator::overlap, atoms);
  //OUT() << S << std::endl;

  LOG(1, "HartreeFock") << "shell: l ao cg" << std::endl;
  for (auto &shell: shells) {
    LOG(1, "HartreeFock") << "shell: "
      << shell.size() << " "
      << shell.nprim() << " "
      << shell.ncontr() << " "
      << std::endl;
  }

  LOG(1, "HartreeFock") << "Calculating kinetic integrals" << std::endl;
  Eigen::MatrixXd T = getOneBodyIntegrals(
    shells,
    libint2::Operator::kinetic,
    atoms);
  LOG(1, "HartreeFock") << "T(" << T.rows() << "," << T.cols() << ")" << std::endl;


  LOG(1, "HartreeFock") << "Compute nuclear repulsion integrals" << std::endl;
  Eigen::MatrixXd V = getOneBodyIntegrals(
    shells,
    libint2::Operator::nuclear,
    atoms);
  LOG(1, "HartreeFock") << "V(" << V.rows() << "," << V.cols() << ")" << std::endl;

  LOG(1, "HartreeFock") << "Calculating the core hamiltonian" << std::endl;
  Eigen::MatrixXd H = T + V;

  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);

  Eigen::MatrixXd D(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd F(nBasisFunctions, nBasisFunctions);
  Eigen::MatrixXd D_last = D;

  D *= 0.0;

  for ( i=0 ; i < nBasisFunctions ; i++) {
    for (unsigned j=i ; j < i+1 ; j++) {
      D(i,j) = 1;
    }
  }

  LOG(1, "HartreeFock") << "Initial Density Matrix:" << std::endl;
  //OUT() << D << std::endl;


  unsigned int iter(0);
  double rmsd(0);
  double energyDifference(0);
  double ehf(0);
  double ehfLast(0);

  do {

    ++iter;

    ehfLast = ehf;
    D_last = D;

    F = H;
    F += getTwoBodyFock(shells, D);

    // solve F C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gen_eig_solver(F, S);
    auto eps(gen_eig_solver.eigenvalues());
    // C now has all eigenvectors from the shell of course
    auto C(gen_eig_solver.eigenvectors());

    // Computer the new density D = C(occ) . C(occ)T
    auto C_occ(C.leftCols(No));
    D = C_occ * C_occ.transpose();

    ehf = 0.0;
    for (unsigned i=0 ; i < nBasisFunctions; i++) {
      for (unsigned j=0 ; j < nBasisFunctions ; j++) {
        ehf += D(i,j) * (H(i,j) + F(i,j));
      }
    }

    energyDifference = ehf - ehfLast;
    rmsd = (D - D_last).norm();

    if (iter == 1) {
      LOG(0, "HartreeFock") <<
        "Iter" << "\t" <<
        "E(elec)" << "\t" <<
        "E(tot)" << "\t" <<
        "Delta(E)" << "\t" <<
        "RMS(D)"
        << std::endl;
    }

    LOG(0, "HartreeFock") << std::setprecision(15) << std::setw(10) <<
        iter << "\t" <<
        ehf  << "\t" <<
        ehf + enuc << "\t" <<
        energyDifference << "\t" <<
        rmsd << "\t" <<
        std::endl;

  } while (
      ((fabs(energyDifference) > electronicConvergence) ||
       (fabs(rmsd) > electronicConvergence))        &&
      (iter < maxIterations)
      );

  LOG(0, "HartreeFock") << "energy=" << ehf + enuc << std::endl;

  libint2::finalize();  // do not use libint after this

}
