#include <string>
#include <vector>
#include <algorithm>
#include <algorithms/CoulombIntegralsFromRotatedCoulombIntegrals.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <ctf.hpp>
#include <numeric>
#include <set>
#include <map>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromRotatedCoulombIntegrals);

struct IntegralProvider {

  IntegralProvider(
    size_t No_,
    size_t Nv_,
    std::vector<double> &coefficients,
    CTF::Tensor<double> &coulombIntegrals)
    :No(No_), Nv(Nv_), Np(No_+Nv_), C(coefficients)
 {
    std::vector<int64_t> indices(Np*Np*Np*Np);
    std::iota(indices.begin(), indices.end(), 0);
    Vklmn.resize(indices.size());
    coulombIntegrals.read(indices.size(), indices.data(), Vklmn.data());
  }

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

  // This is a faster version of the obvious implementation
  std::vector<double> compute_fast(Index P, Index Q, Index R, Index S) {
    // < P Q | R S > = (P R | Q S)

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
  const std::vector<double> &C;
  std::vector<double> Vklmn;
};


CTF::Tensor<double>*
stdVectorToTensor(const std::vector<double> v, const std::vector<int> lens) {
  std::vector<int> syms(lens.size(), NS);
  auto result(
    new CTF::Tensor<double>(
      lens.size(), lens.data(), syms.data(), *Cc4s::world, "T"));
  std::vector<int64_t> indices(v.size());
  std::iota(indices.begin(), indices.end(), 0);
  result->write(v.size(), indices.data(), v.data());
  return result;
}


void CoulombIntegralsFromRotatedCoulombIntegrals::run() {
  auto C(getTensorArgument("OrbitalCoefficients"));
  auto V(getTensorArgument("CoulombIntegrals"));
  const int nelect(getIntegerArgument("nelec", -1));
  const int No(getIntegerArgument("No", nelect/2));
  const int Nv(getIntegerArgument("Nv", nelect/2));
  const int Np(No + Nv);
  std::vector<double> orbitals;

  orbitals.resize(Np*Np);
  {
    std::vector<int64_t> indices(Np*Np);
    std::iota(indices.begin(), indices.end(), 0);
    C->read(orbitals.size(), indices.data(), orbitals.data());
  }

  IntegralProvider engine(No, Nv, orbitals, *V);

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

  LOG(1, "CoulombIntegralsFromRotatedCoulombIntegrals")
    << "No: " << No << std::endl;
  LOG(1, "CoulombIntegralsFromRotatedCoulombIntegrals")
    << "Nv: " << Nv << std::endl;

  for (const auto &integral : integralInfos) {
    if ( ! isArgumentGiven(integral.name) ) continue;
    const auto& i(integral.indices);

    LOG(1, "CoulombIntegralsFromRotatedCoulombIntegrals")
      << "Computing " <<  integral.name << std::endl;

    std::vector<double> result;
    if (!isArgumentGiven("fast")) {
      LOG(1, "CoulombIntegralsFromRotatedCoulombIntegrals")
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
      integral.name, stdVectorToTensor(result, lens));
  }

}
