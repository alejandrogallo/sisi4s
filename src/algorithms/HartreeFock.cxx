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
compute_nuclear_repulsion_energy(std::vector<libint2::Atom>& structure)
{
  int i, j;
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
compute_1body_ints(
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
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

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
compute_2body_fock_simple(
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
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(auto s4=0; s4!=shells.size(); ++s4) {

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
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * buf_1234[f1234];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);
          const auto* buf_1324 = buf[0];

          for(auto f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1324) {
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
  LOG(0, "HartreeFock") << "Calculating hartree fock" << std::endl;
  LOG(1, "HartreeFock") << "We'll see how this works" << std::endl;

  const std::string xyzStructureFile(
    getTextArgument("xyzStructureFile", "")
  );
  const std::string basis_set(
    getTextArgument("basisSet", "sto-3g")
  );
  int electronic_iterations(
    getIntegerArgument("maxIterations", 16)
  );
  double electronic_convergence(
    getRealArgument("energyDifference", 1e-4)
  );
  int number_of_electrons(
    getIntegerArgument("numberOfElectrons", 0)
  );
  int i;
  int n_basis_functions;
  int No;
  double enuc;

  LOG(1, "HartreeFock") << "Electronic iterations  = " << electronic_iterations  << std::endl;
  LOG(1, "HartreeFock") << "Electronic convergence = " << electronic_convergence << std::endl;
  LOG(1, "HartreeFock") << "Basis set              = " << basis_set              << std::endl;

  LOG(2, "HartreeFock") << "Initialising libint2" << std::endl;
  libint2::initialize();  // safe to use libint now

  LOG(1, "HartreeFock") << "Reading structure from " << xyzStructureFile << std::endl;
  std::ifstream structureFileStream(xyzStructureFile.c_str());
  std::vector<libint2::Atom> atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();

  LOG(1, "HartreeFock") << "Computing number of electrons" << std::endl;
  for (i = 0; i < atoms.size(); ++i) {
    number_of_electrons += atoms[i].atomic_number;
  }
  LOG(1, "HartreeFock") << "Number of atoms = " << atoms.size() << std::endl;
  LOG(1, "HartreeFock") << "Number of electrons = " << number_of_electrons << std::endl;

  // restricted hartree fock
  No = number_of_electrons/2;

  enuc = compute_nuclear_repulsion_energy(atoms);
  LOG(1, "HartreeFock") << "Nuclear repulsion energy " << enuc << std::endl;

  LOG(1, "HartreeFock") << "Initializing basis set" << std::endl;
  libint2::BasisSet shells(basis_set, atoms);

  n_basis_functions = shells.nbf();

  LOG(1, "HartreeFock") << "Number of basis functions = " << n_basis_functions
    << std::endl;
  LOG(1, "HartreeFock")
    << "Maximum number of primitives in shells = "
    << shells.max_nprim()
    << std::endl;
  LOG(1, "HartreeFock")
    << "Maximum value of l in shells = "
    << shells.max_l()
    << std::endl;

  LOG(1, "HartreeFock") << "Calculating overlaps" << std::endl;
  Eigen::MatrixXd S = compute_1body_ints(
    shells, libint2::Operator::overlap, atoms
  );
  OUT() << S << std::endl;

  LOG(1, "HartreeFock") << "Calculating kinteic integrals" << std::endl;
  Eigen::MatrixXd T = compute_1body_ints(
    shells, libint2::Operator::kinetic, atoms
  );
  OUT() << T << std::endl;

  LOG(1, "HartreeFock") << "Compute nuclear repulsion integrals" << std::endl;
  Eigen::MatrixXd V = compute_1body_ints(
    shells, libint2::Operator::nuclear, atoms
  );
  OUT() << V << std::endl;

  LOG(1, "HartreeFock") << "Calculating the core hamiltonian" << std::endl
    << "  H = T + V" << std::endl;
  Eigen::MatrixXd H = T + V;

  // T and V no longer needed, free up the memory
  T.resize(0,0);
  V.resize(0,0);

  Eigen::MatrixXd D(n_basis_functions, n_basis_functions);
  Eigen::MatrixXd F(n_basis_functions, n_basis_functions);
  Eigen::MatrixXd D_last = D;

  D *= 0.0;

  for ( i=0 ; i < n_basis_functions ; i++) {
    for (unsigned j=i ; j < i+1 ; j++) {
      D(i,j) = 1;
    }
  }

  LOG(1, "HartreeFock") << "Initial Density Matrix:" << std::endl;
  OUT() << D << std::endl;


  int iter(0);
  double rmsd(0);
  double energy_diff(0);
  double ehf(0);
  double ehf_last(0);

  do {

    ++iter;

    ehf_last = ehf;
    D_last = D;

    F = H;
    F += compute_2body_fock_simple(shells, D);

    // solve F C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gen_eig_solver(F, S);
    auto eps(gen_eig_solver.eigenvalues());
    // C now has all eigenvectors from the shell of course
    auto C(gen_eig_solver.eigenvectors());

    // Computer the new density D = C(occ) . C(occ)T
    auto C_occ(C.leftCols(No));
    D = C_occ * C_occ.transpose();

    ehf = 0.0;
    for (unsigned i=0 ; i < n_basis_functions; i++) {
      for (unsigned j=0 ; j < n_basis_functions ; j++) {
        ehf += D(i,j) * (H(i,j) + F(i,j));
      }
    }

    energy_diff = ehf - ehf_last;
    rmsd = (D - D_last).norm();

    if (iter == 1)
      LOG(0, "HartreeFock") <<
        "Iter   "
        "E(elec)   "
        "E(tot)    "
        "Delta(E)   "
        "RMS(D)  "
        << std::endl;

    LOG(0, "HartreeFock") <<
        iter << "      " <<
        ehf  << "   " <<
        ehf + enuc << "    " <<
        energy_diff << "    " <<
        rmsd << "    " <<
        std::endl;

  } while (
      ((fabs(energy_diff) > electronic_convergence) ||
       (fabs(rmsd) > electronic_convergence))        &&
      (iter < electronic_iterations)
      );

  LOG(0, "HartreeFock") << "energy=" << ehf + enuc << std::endl;

  libint2::finalize();  // do not use libint after this

}
