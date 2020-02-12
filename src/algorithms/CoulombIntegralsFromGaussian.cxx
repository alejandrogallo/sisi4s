#include <string>
#include <vector>
#include <algorithm>
#include <libint2.hpp>
#include <algorithms/CoulombIntegralsFromGaussian.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromGaussian);

// struct for storing information about the shell ends in the for loops
// calculating the integrals
struct ShellInfo {
  // size of the shell, global begin and global end
  // l should be the angular momentum
  size_t size, begin, end, l;
  // constructor from a BasisSet and a shell index i
  inline ShellInfo(const libint2::BasisSet &shells, const size_t i) {
    size = shells[i].size();
    begin = shells.shell2bf()[i];
    end = size + begin;
    l = shells[i].contr[0].l;
  }
};

struct CoulombIntegralsProvider {
  CoulombIntegralsProvider(
    const size_t No_,
    const size_t Nv_,
    const libint2::BasisSet& shells_):
    No(No_), Nv(Nv_), Np(No_+Nv_), shells(shells_) {}

  void compute() {
    // If already computed, return
    libint2::initialize();

    const size_t NpNpNpNp(Np*Np*Np*Np);
    LOG(1, "Integrals")
      << "Allocating and computing Vklmn ("
      << sizeof(double) * NpNpNpNp / std::pow(2, 30)
      << " GB)" << std::endl;
    Vklmn.resize(NpNpNpNp, 0.0);

    libint2::Engine engine(
      libint2::Operator::coulomb,
      shells.max_nprim(),
      shells.max_l(), 0);

    // store shell by shell calculation in this buffer
    const auto& vsrqp = engine.results();

    // the outside loops will loop over the shells.
    // This will create a block of Vpqrs, where pqrs are contracted
    // gaussian indices belonging to their respective shells.
    for (size_t _K(0); _K < shells.size(); ++_K) { // kappa
    for (size_t _L(0); _L < shells.size(); ++_L) { // lambda
    for (size_t _M(0); _M < shells.size(); ++_M) { // mu
    for (size_t _N(0); _N < shells.size(); ++_N) { // nu
      const ShellInfo K(shells, _K),
                      L(shells, _L),
                      M(shells, _M),
                      N(shells, _N);

      // compute integrals (K L , M N)
      engine.compute(shells[_K], shells[_L], shells[_M], shells[_N]);

        for (size_t k(K.begin), Inmlk = 0; k < K.end; ++k) {
        for (size_t l(L.begin); l < L.end; ++l) {
        for (size_t m(M.begin); m < M.end; ++m) {
        for (size_t n(N.begin); n < N.end; ++n, ++Inmlk) {

          // <p q | r s> = (p r , q s)

          size_t bigI(
            n +
            m * Np +
            l * Np*Np +
            k * Np*Np*Np);

          Vklmn[bigI] += vsrqp[0][Inmlk];

        } // n
        } // m
        } // l
        } // k

    } // N
    } // M
    } // L
    } // K

  libint2::finalize();

  }

  double *data() {return Vklmn.data();}

  private:
  const size_t No, Nv, Np;
  const libint2::BasisSet& shells;
  std::vector<double> Vklmn;
};


void CoulombIntegralsFromGaussian::run() {
  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""));
  const std::string basisSet(getTextArgument("basisSet"));
  const int nelect(getIntegerArgument("nelec", -1));
  const int No(getIntegerArgument("No", nelect/2));

  std::ifstream structureFileStream(xyzStructureFile.c_str());
  const auto atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();
  const libint2::BasisSet shells(basisSet, atoms);
  const int Np(shells.nbf());
  const int Nv(Np - No);

  CoulombIntegralsProvider engine(No, Nv, shells);

  LOG(1, "CoulombIntegralsFromGaussian")
    << "structure: " << xyzStructureFile << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "No: " << No << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "Nv: " << Nv << std::endl;

  const std::vector<int> lens(4, Np);
  const std::vector<int> syms(4, NS);
  auto Vklmn(new CTF::Tensor<double>(4, lens.data(), syms.data(),
                                      *Cc4s::world, "V"));
  engine.compute();
  {
    std::vector<int64_t> indices(Np*Np*Np*Np);
    std::iota(indices.begin(), indices.end(), 0);
    Vklmn->write(indices.size(), indices.data(), engine.data());
  }
  LOG(1, "CoulombIntegralsFromGaussian")
    << "Allocating CoulombIntegrals" << std::endl;
  allocatedTensorArgument<double>("CoulombIntegrals", Vklmn);

}
