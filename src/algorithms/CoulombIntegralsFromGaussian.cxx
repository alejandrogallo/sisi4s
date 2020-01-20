#include <string>
#include <vector>
#include <libint2.hpp>
#include <algorithms/CoulombIntegralsFromGaussian.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <ctf.hpp>
#include <numeric>      // std::iota

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


struct IntegralProvider {
  IntegralProvider(
    size_t No_,
    size_t Nv_,
    double *coefficients,
    libint2::BasisSet& shells_):
    No(No_), Nv(Nv_), Np(No_+Nv_), C(coefficients), shells(shells_) {}

  enum Index { NO, NV, NP };
  struct Limit {
    size_t lower; size_t upper; size_t size;
    Limit(size_t l, size_t h): lower(l), upper(h), size(h - l){}
  };
  Limit indexToLimits(Index &i) {
    if      (i == Index::NO) return Limit(0, No);
    else if (i == Index::NV) return Limit(No, Nv);
    else                     return Limit(0, Nv);
  }

  void compute_Vklmn() {
    // If already computed, return
    if (Vklmn) return;

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

          // TODO: do it correctly
          Vklmn[Inmlk] += vsrqp[0][Inmlk];

        } // s
        } // r
        } // q
        } // p

    } // N
    } // M
    } // L
    } // K

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
    for (size_t p(qLim.lower), Isrqp(0); p < qLim.upper; ++p         ) {
    for (size_t q(qLim.lower);           q < qLim.upper; ++q         ) {
    for (size_t r(rLim.lower);           r < rLim.upper; ++r         ) {
    for (size_t s(sLim.lower);           s < sLim.upper; ++s, ++Isrqp) {

      for (size_t k(0), Inmlk = 0; k < Np; ++k         ) {
      for (size_t l(0);            l < Np; ++l         ) {
      for (size_t m(0);            m < Np; ++m         ) {
      for (size_t n(0);            n < Np; ++n, ++Inmlk) {

        // <p q | r s> = (p r , q s)
        //               (k l , m n)

        result.data()[Isrqp] +=
          C[k + p*Np] *
          C[l + r*Np] *
          C[m + q*Np] *
          C[n + s*Np] *
          Vklmn[Inmlk];

      } // s
      } // r
      } // q
      } // p

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


void CoulombIntegralsFromGaussian::run() {
  auto C(getTensorArgument("OrbitalCoefficients"));
  int antisymmetrize(getIntegerArgument("antisymmetrize", 0));
  isArgumentGiven("PHPHCoulombIntegrals");
  isArgumentGiven("PPHHCoulombIntegrals");
  isArgumentGiven("HHHHCoulombIntegrals");
  isArgumentGiven("HHHPCoulombIntegrals");
  isArgumentGiven("PPPPCoulombIntegrals");
  isArgumentGiven("PPPHCoulombIntegrals");
  isArgumentGiven("PHHHCoulombIntegrals");
  isArgumentGiven("HHPPCoulombIntegrals");
  isArgumentGiven("PHHPCoulombIntegrals");
  isArgumentGiven("HPHHCoulombIntegrals");
  isArgumentGiven("HPHPCoulombIntegrals");
  isArgumentGiven("HPPPCoulombIntegrals");
  isArgumentGiven("PPHPCoulombIntegrals");
  isArgumentGiven("HPPHCoulombIntegrals");
  isArgumentGiven("HHPHCoulombIntegrals");
  isArgumentGiven("PHPPCoulombIntegrals");
}
