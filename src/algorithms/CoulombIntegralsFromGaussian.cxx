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
    inline const size_t operator[](const size_t& i) const { return i - lower; }
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

  // This is a faster version of the obvious implementation
  std::vector<double> compute_fast(Index P, Index Q, Index R, Index S) {
    // < P Q | R S > = (P R | Q S)

    compute_Vklmn();

    const Limit pLim(indexToLimits(P)),
                qLim(indexToLimits(Q)),
                rLim(indexToLimits(R)),
                sLim(indexToLimits(S));

    const size_t dimension(pLim.size * qLim.size * rLim.size * sLim.size);
    std::vector<double> result(dimension, 0.0);

    for (size_t s(sLim.lower); s < sLim.upper; ++s) {
    for (size_t q(qLim.lower); q < qLim.upper; ++q) {
    for (size_t r(rLim.lower); r < rLim.upper; ++r) {
    for (size_t p(pLim.lower); p < pLim.upper; ++p) {

      //  <p q | r s> = (p r , q s)
      const int IPQRS(
        pLim[p] +
        qLim[q] * pLim.size +
        rLim[r] * pLim.size*qLim.size +
        sLim[s] * pLim.size*qLim.size*rLim.size);

      std::vector<double> Vpmlk(
        Np*Np*Np* pLim.size, 0.0);
      std::vector<double> Vprlk(
        Np*Np*    rLim.size*pLim.size, 0.0);
      std::vector<double> Vprqk(
        Np*       qLim.size*rLim.size*pLim.size, 0.0);

      // we will define four intermediates
      for (size_t k(0), Inmlk = 0; k < Np; ++k) {

        const size_t Iprqk(
          pLim[p] +
          rLim[r] * pLim.size +
          qLim[q] * pLim.size*rLim.size +
          k       * pLim.size*rLim.size*pLim.size);

        result.data()[IPQRS] += C[k + s*Np] * Vprqk[Iprqk];

      for (size_t l(0); l < Np; ++l) {

        const size_t Iprlk(
          pLim[p] +
          rLim[r] * pLim.size +
          l       * pLim.size*rLim.size +
          k       * pLim.size*rLim.size*Np);

        Vprqk[Iprqk] += C[l + q*Np] * Vprlk[Iprlk];

      for (size_t m(0); m < Np; ++m) {

        const size_t Ipmlk(
          pLim[p] +
          m * pLim.size +
          l * pLim.size*Np +
          k * pLim.size*Np*Np);

        Vprlk[Iprlk] += C[m + r*Np] * Vpmlk[Ipmlk];

      for (size_t n(0); n < Np; ++n, ++Inmlk) {

        //  Vnmlk
        Vpmlk[Ipmlk] += C[n + p*Np] * Vklmn[Inmlk];

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

  std::vector<double> compute(Index P, Index Q, Index R, Index S) {
    // < P Q | R S > = (P R | Q S)

    compute_Vklmn();

    const Limit pLim(indexToLimits(P)),
                qLim(indexToLimits(Q)),
                rLim(indexToLimits(R)),
                sLim(indexToLimits(S));

    const size_t dimension(pLim.size * qLim.size * rLim.size * sLim.size);
    std::vector<double> result(dimension, 0.0);

    for (size_t s(sLim.lower), ipqrs(0); s < sLim.upper; ++s         ) {
    for (size_t q(qLim.lower);           q < qLim.upper; ++q         ) {
    for (size_t r(rLim.lower);           r < rLim.upper; ++r         ) {
    for (size_t p(pLim.lower);           p < pLim.upper; ++p, ++ipqrs) {

      for (size_t k(0), Inmlk = 0; k < Np; ++k         ) {
      for (size_t l(0);            l < Np; ++l         ) {
      for (size_t m(0);            m < Np; ++m         ) {
      for (size_t n(0);            n < Np; ++n, ++Inmlk) {

        // <p q | r s> = (p r , q s)
        //LOG(1, "Integrals:")  << ipqrs  << ": " << Inmlk << std::endl;
        const int IPQRS(
          pLim[p] +
          qLim[q] * pLim.size +
          rLim[r] * pLim.size*qLim.size +
          sLim[s] * pLim.size*qLim.size*rLim.size);

        result.data()[IPQRS] +=
          C[k + s*Np] *
          C[l + q*Np] *
          C[m + r*Np] *
          C[n + p*Np] *
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
  const int nelect(getIntegerArgument("nelec", -1));
  const int No(getIntegerArgument("No", nelect/2));
  int Nv, Np;
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


  IntegralProvider engine(No, Nv, orbitals.data(), shells);

  const std::vector<IntegralInfo> integralInfos({
    {"HHHHCoulombIntegrals", {NO,NO,NO,NO}, "ijkl"},
    {"HHHPCoulombIntegrals", {NO,NO,NO,NV}, "ijka"},
    {"HHPHCoulombIntegrals", {NO,NO,NV,NO}, "ijak"},
    {"HHPPCoulombIntegrals", {NO,NO,NV,NV}, "ijab"},
    {"HPHHCoulombIntegrals", {NO,NV,NO,NO}, "iajk"},
    {"HPHPCoulombIntegrals", {NO,NV,NO,NV}, "iajb"},
    {"HPPHCoulombIntegrals", {NO,NV,NV,NO}, "iabj"},
    {"HPPPCoulombIntegrals", {NO,NV,NV,NV}, "iabc"},
    {"PHHHCoulombIntegrals", {NV,NO,NO,NO}, "aijk"},
    {"PHHPCoulombIntegrals", {NV,NO,NO,NV}, "aijb"},
    {"PHPHCoulombIntegrals", {NV,NO,NV,NO}, "aibj"},
    {"PHPPCoulombIntegrals", {NV,NO,NV,NV}, "aibc"},
    {"PPHHCoulombIntegrals", {NV,NV,NO,NO}, "abij"},
    {"PPHPCoulombIntegrals", {NV,NV,NO,NV}, "abic"},
    {"PPPHCoulombIntegrals", {NV,NV,NV,NO}, "abci"},
    {"PPPPCoulombIntegrals", {NV,NV,NV,NV}, "abcd"},
  });

  LOG(1, "CoulombIntegralsFromGaussian")
    << "structure: " << xyzStructureFile << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "No: " << No << std::endl;
  LOG(1, "CoulombIntegralsFromGaussian") << "Nv: " << Nv << std::endl;

  for (const auto &integral : integralInfos) {
    if ( ! isArgumentGiven(integral.name) ) continue;
    const auto& i(integral.indices);

    LOG(1, "CoulombIntegralsFromGaussian")
      << "Computing " <<  integral.name << std::endl;

    std::vector<double> result;
    if (!isArgumentGiven("fast")) {
      LOG(1, "CoulombIntegralsFromGaussian")
        << "Computing slow implementation" << std::endl;
      result = engine.compute(i[0], i[1], i[2], i[3]);
    } else {
      result = engine.compute_fast(i[0], i[1], i[2], i[3]);
    }
    std::vector<int> lens(4);
    for (unsigned int j(0) ; j < 4 ; j++) {
      lens[j] = (int)engine.indexToInt(i[j]);
    }
    allocatedTensorArgument<double>(
      integral.name, vectorToTensor(result, lens));
  }

}
