#include <string>
#include <vector>
#include <algorithm>
#include <libint2.hpp>
#include <algorithms/CoulombIntegralsFromGaussian.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <ctf.hpp>
#include <numeric>      // std::iota
#include <map>      // std::iota

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromGaussian);

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

enum Index { NO, NV, NP };
typedef std::array<Index, 4> PQRSIndices;

struct IntegralProvider {
  IntegralProvider(
    size_t No_,
    size_t Nv_,
    double *coefficients,
    libint2::BasisSet& shells_):
    No(No_), Nv(Nv_), Np(No_+Nv_), C(coefficients), shells(shells_) {}

  struct Limit {
    size_t lower; size_t upper; size_t size;
    Limit(size_t l, size_t h): lower(l), upper(h), size(h - l){}
  };
  Limit indexToLimits(Index &i) {
    if      (i == Index::NO) return Limit(0 , No     );
    else if (i == Index::NV) return Limit(No, No + Nv);
    else                     return Limit(0 , No + Nv);
  }
  size_t indexToInt(Index i) {
    if      (i == Index::NO) return No;
    else if (i == Index::NV) return Nv;
    else                     return No + Nv;
  }

  void compute_Vklmn() {
    // If already computed, return
    if (Vklmn) return;
    libint2::initialize();

    const size_t NpNpNpNp(Np*Np*Np*Np);
    LOG(1, "Integrals")
      << "Allocating and computing Vklmn ("
      << sizeof(double) * NpNpNpNp / std::pow(2, 30)
      << " GB)" << std::endl;
    Vklmn = new double[NpNpNpNp];

    libint2::Engine engine(
      libint2::Operator::coulomb,
      shells.max_nprim(),
      shells.max_l(), 0);

    // store shell by shell calculation in this buffer
    const auto& vsrqp = engine.results();

    for (size_t p(0); p < NpNpNpNp; ++p) { Vklmn[p] = 0.0; }

    LOG(1, "Integrals") << "Initialized to zero" << std::endl;

    // the outside loops will loop over the shells.
    // This will create a block of Vpqrs, where pqrs are contracted
    // gaussian indices belonging to their respective shells.
    // Since we only want to calculate integrals transformed by the
    // coefficients provided, we will contract them with the coefficients.
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
          //               (k l , m n)

          // TODO: do it easier
          size_t bigI(
            n +
            m * Np +
            l * Np*Np +
            k * Np*Np*Np);

          //LOG(1, "Integrals")  << Inmlk << " <- "<< bigI << std::endl;

          Vklmn[bigI] += vsrqp[0][Inmlk];

        } // s
        } // r
        } // q
        } // p

    } // N
    } // M
    } // L
    } // K

  libint2::finalize();

  }

  std::vector<double> compute(Index P, Index Q, Index R, Index S) {

    compute_Vklmn();

    const Limit pLim(indexToLimits(P)),
                qLim(indexToLimits(Q)),
                rLim(indexToLimits(R)),
                sLim(indexToLimits(S));

    const size_t dimension(pLim.size * qLim.size * rLim.size * sLim.size);
    std::vector<double> result(dimension, 0.0);

    // loop over all C orbitals
    for (size_t p(pLim.lower), Isrqp(0); p < pLim.upper; ++p         ) {
    for (size_t r(rLim.lower);           r < rLim.upper; ++r         ) {
    for (size_t q(qLim.lower);           q < qLim.upper; ++q         ) {
    for (size_t s(sLim.lower);           s < sLim.upper; ++s, ++Isrqp) {

      for (size_t k(0), Inmlk = 0; k < Np; ++k         ) {
      for (size_t l(0);            l < Np; ++l         ) {
      for (size_t m(0);            m < Np; ++m         ) {
      for (size_t n(0);            n < Np; ++n, ++Inmlk) {

        // <p q | r s> = (p r , q s)
        //LOG(1, "Integrals:")  << Isrqp  << ": " << Inmlk << std::endl;

        result.data()[Isrqp] +=
          C[k + p*Np] *
          C[l + r*Np] *
          C[m + q*Np] *
          C[n + s*Np] *
          Vklmn[Inmlk];

      } // n
      } // m
      } // l
      } // k

    } // s
    } // r
    } // q
    } // p

    return result;

  }

  private:
  size_t No, Nv, Np;
  double *C;
  libint2::BasisSet& shells;
  double *Vklmn = nullptr;
};


CTF::Tensor<double>*
vectorToTensor(const std::vector<double> v, const std::vector<int> lens) {
  std::vector<int> syms(lens.size(), NS);
  auto result(
    new CTF::Tensor<double>(
      lens.size(), lens.data(), syms.data(), *Cc4s::world, "T"));
  std::vector<int64_t> indices(v.size());
  std::iota(indices.begin(), indices.end(), 0);
  result->write(v.size(), indices.data(), v.data());
  return result;
}


void CoulombIntegralsFromGaussian::run() {
  const std::string xyzStructureFile(getTextArgument("xyzStructureFile", ""));
  const std::string basisSet(getTextArgument("basisSet"));
  auto C(getTensorArgument("OrbitalCoefficients"));
  int No(getIntegerArgument("No")), Nv, Np;
  std::vector<double> orbitals;

  std::ifstream structureFileStream(xyzStructureFile.c_str());
  const auto atoms(libint2::read_dotxyz(structureFileStream));
  structureFileStream.close();
  libint2::BasisSet shells(basisSet, atoms);

  Np = shells.nbf();
  Nv = Np - No;

  orbitals.resize(Np*Np);
  {
    std::vector<int64_t> indices(Np*Np);
    std::iota(indices.begin(), indices.end(), 0);
    C->read(orbitals.size(), indices.data(), orbitals.data());
  }

  IntegralProvider ints(No, Nv, orbitals.data(), shells);

  std::map< std::string, std::array<Index, 4> > integralSizes({
      { "HHHHCoulombIntegrals", PQRSIndices({NO, NO, NO, NO})},
      { "HHHPCoulombIntegrals", PQRSIndices({NO, NO, NO, NV})},
      { "HHPHCoulombIntegrals", PQRSIndices({NO, NO, NV, NO})},
      { "HHPPCoulombIntegrals", PQRSIndices({NO, NO, NV, NV})},
      { "HPHHCoulombIntegrals", PQRSIndices({NO, NV, NO, NO})},
      { "HPHPCoulombIntegrals", PQRSIndices({NO, NV, NO, NV})},
      { "HPPHCoulombIntegrals", PQRSIndices({NO, NV, NV, NO})},
      { "HPPPCoulombIntegrals", PQRSIndices({NO, NV, NV, NV})},
      { "PHHHCoulombIntegrals", PQRSIndices({NV, NO, NO, NO})},
      { "PHHPCoulombIntegrals", PQRSIndices({NV, NO, NO, NV})},
      { "PHPHCoulombIntegrals", PQRSIndices({NV, NO, NV, NO})},
      { "PHPPCoulombIntegrals", PQRSIndices({NV, NO, NV, NV})},
      { "PPHHCoulombIntegrals", PQRSIndices({NV, NV, NO, NO})},
      { "PPHPCoulombIntegrals", PQRSIndices({NV, NV, NO, NV})},
      { "PPPHCoulombIntegrals", PQRSIndices({NV, NV, NV, NO})},
      { "PPPPCoulombIntegrals", PQRSIndices({NV, NV, NV, NV})},
  });

  LOG(1, "CoulombIntegralsFromGaussian")
    << "structure: " << xyzStructureFile << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "No: " << No << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "Nv: " << Nv << std::endl;

  for (const auto &integral : integralSizes) {
    const auto& name(integral.first);
    if (isArgumentGiven(name)) {
      const auto& i(integral.second);

      LOG(1, "CoulombIntegralsFromGaussian")
        << "Computing " <<  name << std::endl;

      std::vector<double> result(ints.compute(i[0], i[1], i[2], i[3]));
      std::vector<int> lens(4);
      for (unsigned int j(0) ; j < 4 ; j++) {
        lens[j] = (int)ints.indexToInt(i[j]);
      }
      allocatedTensorArgument<double>(name, vectorToTensor(result, lens));
    }
  }

}
