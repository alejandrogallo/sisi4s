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
#include <util/Emitter.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromRotatedCoulombIntegrals);
#define LOGGER(_l) LOG(_l, "CoulombIntegralsFromRotatedCoulombIntegrals")

template <typename V>
struct IntegralProvider {

  IntegralProvider(size_t no, size_t nv, bool chemistNotation_)
    : No(no), Nv(nv), Np(no+nv), chemistNotation(chemistNotation_) {}

  struct Limit {
    size_t lower, upper, size;
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

  virtual V compute(Index P, Index Q, Index R, Index S) = 0;

  const size_t No, Nv, Np;
  const bool chemistNotation = false;
};

struct VectorIntegralProvider: public IntegralProvider< std::vector<double> > {

  VectorIntegralProvider( size_t no
                        , size_t nv
                        , bool chemistNotation_
                        , CTF::Tensor<double> &coeffs
                        , CTF::Tensor<double> &coulombIntegrals
                        ) :IntegralProvider(no, nv, chemistNotation_)
  {
    if (! chemistNotation_ )
      throw EXCEPTION("Physics notation not supported for this provider");

    const int rank_m = int(Cc4s::world->rank == 0); // rank mask

    // read in the orbital coefficients
    C.resize(rank_m * Np*Np);
    std::vector<int64_t> indices(C.size());
    std::iota(indices.begin(), indices.end(), 0);
    coeffs.read(C.size(), indices.data(), C.data());

    // read in the integrals
    indices.resize(rank_m * Np*Np*Np*Np);
    std::iota(indices.begin(), indices.end(), 0);
    Vklmn.resize(indices.size());
    coulombIntegrals.read(indices.size(), indices.data(), Vklmn.data());
  }

  std::vector<double> compute(Index P, Index Q, Index R, Index S) {

    const Limit pLim(indexToLimits(P))
              , qLim(indexToLimits(Q))
              , rLim(indexToLimits(R))
              , sLim(indexToLimits(S))
              ;

    const size_t dimension(pLim.size * qLim.size * rLim.size * sLim.size);
    std::vector<double> result(dimension, 0.0);

    std::vector<double> Vpmlk(Np*Np*Np* pLim.size, 0.0)
                      , Vprlk(Np*Np*    pLim.size*rLim.size, 0.0)
                      , Vprqk(Np*       pLim.size*rLim.size*qLim.size, 0.0)
                      ;

    // COMPUTE Vpmlk =====================
    for (size_t p(pLim.lower); p < pLim.upper; ++p) {
      for (size_t k(0), Inmlk = 0; k < Np; ++k) {
      for (size_t l(0); l < Np; ++l) {
      for (size_t m(0); m < Np; ++m) {

        const size_t Ipmlk = pLim[p]
                           + m * pLim.size
                           + l * pLim.size*Np
                           + k * pLim.size*Np*Np
                           ;

      for (size_t n(0); n < Np; ++n, ++Inmlk) {

        Vpmlk[Ipmlk] += C[n + p*Np] * Vklmn[Inmlk];

      } /* n */ } /* m */ } /* l */ } /* k */
    }

    // COMPUTE Vprlk ======================
    for (size_t r(rLim.lower); r < rLim.upper; ++r) {
    for (size_t p(pLim.lower); p < pLim.upper; ++p) {
      for (size_t k(0); k < Np; ++k) {
      for (size_t l(0); l < Np; ++l) {

        const size_t Iprlk = pLim[p]
                           + rLim[r] * pLim.size
                           + l       * pLim.size*rLim.size
                           + k       * pLim.size*rLim.size*Np
                           ;


      for (size_t m(0); m < Np; ++m) {

        const size_t Ipmlk = pLim[p]
                           + m * pLim.size
                           + l * pLim.size*Np
                           + k * pLim.size*Np*Np
                           ;

        Vprlk[Iprlk] += C[m + r*Np] * Vpmlk[Ipmlk];

      } /* m */ } /* l */ } /* k */
    }}

    for (size_t r(rLim.lower); r < rLim.upper; ++r) {
    for (size_t q(qLim.lower); q < qLim.upper; ++q) {
    for (size_t p(pLim.lower); p < pLim.upper; ++p) {
      // COMPUTE Vprqk =====================
      for (size_t k(0); k < Np; ++k) {

        const size_t Iprqk = pLim[p]
                           + rLim[r] * pLim.size
                           + qLim[q] * pLim.size*rLim.size
                           + k       * pLim.size*rLim.size*qLim.size
                           ;

      for (size_t l(0); l < Np; ++l) {

        const size_t Iprlk = pLim[p] +
                             rLim[r] * pLim.size +
                             l       * pLim.size*rLim.size +
                             k       * pLim.size*rLim.size*Np
                             ;

        Vprqk[Iprqk] += C[l + q*Np] * Vprlk[Iprlk];

      }}
    }}}

    // COMPUTE Vpqrs =====================
    for (size_t s(sLim.lower); s < sLim.upper; ++s) {
    for (size_t r(rLim.lower); r < rLim.upper; ++r) {
    for (size_t q(qLim.lower); q < qLim.upper; ++q) {
    for (size_t p(pLim.lower); p < pLim.upper; ++p) {

      const size_t IPQRS = pLim[p]
                         + qLim[q] * pLim.size
                         + rLim[r] * pLim.size*qLim.size
                         + sLim[s] * pLim.size*qLim.size*rLim.size
                         ;

      for (size_t k(0); k < Np; ++k) {

        const size_t Iprqk = pLim[p]
                           + rLim[r] * pLim.size
                           + qLim[q] * pLim.size*rLim.size
                           + k       * pLim.size*rLim.size*qLim.size
                           ;

        result[IPQRS] += C[k + s*Np] * Vprqk[Iprqk];

      } // k

    }}}}

    return result;

  }
  std::vector<double> C;
  std::vector<double> Vklmn;
};

struct SlowVectorIntegralProvider: public VectorIntegralProvider {
  using VectorIntegralProvider::VectorIntegralProvider;
  std::vector<double> compute(Index P, Index Q, Index R, Index S) {
    // < P Q | R S > = (P R | Q S)

    const Limit pLim(indexToLimits(P))
              , qLim(indexToLimits(Q))
              , rLim(indexToLimits(R))
              , sLim(indexToLimits(S))
              ;

    const size_t dimension(pLim.size * qLim.size * rLim.size * sLim.size);
    std::vector<double> result(dimension, 0.0);

    for (size_t s(sLim.lower), ipqrs(0); s < sLim.upper; ++s         ) {
    for (size_t r(rLim.lower);           r < rLim.upper; ++r         ) {
    for (size_t q(qLim.lower);           q < qLim.upper; ++q         ) {
    for (size_t p(pLim.lower);           p < pLim.upper; ++p, ++ipqrs) {

      for (size_t k(0), Inmlk = 0; k < Np; ++k         ) {
      for (size_t l(0);            l < Np; ++l         ) {
      for (size_t m(0);            m < Np; ++m         ) {
      for (size_t n(0);            n < Np; ++n, ++Inmlk) {

        // <p q | r s> = (p r , q s)
        const int IPQRS = pLim[p]
                        + qLim[q] * pLim.size
                        + rLim[r] * pLim.size * qLim.size
                        + sLim[s] * pLim.size * qLim.size * rLim.size
                        ;

        result[IPQRS] += C[k + s*Np]
                       * C[l + q*Np]
                       * C[m + r*Np]
                       * C[n + p*Np]
                       * Vklmn[Inmlk]
                       ;

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

};

struct CtfIntegralProvider: public IntegralProvider< CTF::Tensor<double> > {
  CtfIntegralProvider(size_t no
                    , size_t nv
                    , bool chemistNotation_
                    , CTF::Tensor<double> &coefficients
                    , CTF::Tensor<double> &coulombIntegrals
                    ):IntegralProvider(no, nv, chemistNotation_)
                    , C(coefficients) , V(coulombIntegrals) {}

  CTF::Tensor<double> *compute() {
    if (VTransformed != nullptr) { return VTransformed; }
    LOGGER(1) << "computing main transformation" << std::endl;
    VTransformed = new CTF::Tensor<double>(true, V);
    if (chemistNotation) {
      (*VTransformed)["pqrs"] = C["ns"]
                              * C["lr"]
                              * C["mq"]
                              * C["kp"]
                              * V["klmn"]
                              ;
    } else {
      (*VTransformed)["pqrs"] = C["ns"]
                              * C["mr"]
                              * C["lq"]
                              * C["kp"]
                              * V["klmn"]
                              ;
    }
    return VTransformed;
  }

  CTF::Tensor<double> compute(Index P, Index Q, Index R, Index S) {
    // Compute the transformation if it's not done
    const auto& vpqrs(compute());

    const Limit p(indexToLimits(P))
              , q(indexToLimits(Q))
              , r(indexToLimits(R))
              , s(indexToLimits(S))
              ;

    int a[] = {(int)p.lower, (int)q.lower, (int)r.lower, (int)s.lower}
      , e[] = {(int)p.upper, (int)q.upper, (int)r.upper, (int)s.upper}
      ;

    LOGGER(1) << "Slicing "
     "{" << a[0] << "," << a[1] << "," << a[2] << "," << a[3] << "}" <<
     " -> " <<
     "{" << e[0] << "," << e[1] << "," << e[2] << "," << e[3] << "}" <<
     std::endl;

    return vpqrs->slice(a, e);
  }

  private:
    CTF::Tensor<double> &C, &V;
    CTF::Tensor<double>* VTransformed = nullptr;
};

CTF::Tensor<double>* stdVectorToTensor( const std::vector<double> v
                                      , const std::vector<int> lens) {
  std::vector<int> syms(lens.size(), NS);
  auto result(new CTF::Tensor<double>( lens.size()
                                     , lens.data()
                                     , syms.data()
                                     , *Cc4s::world
                                     , "T"
                                     ));
  std::vector<int64_t> indices(int(Cc4s::world->rank == 0) * v.size());
  std::iota(indices.begin(), indices.end(), 0);
  result->write(indices.size(), indices.data(), v.data());
  return result;
}

template <typename V>
void computeAndExport ( Algorithm&
                      , IntegralProvider<V>&
                      , std::vector<IntegralInfo>&
                      );

template<>
void computeAndExport ( Algorithm &a
                      , IntegralProvider< std::vector<double> > &engine
                      , std::vector<IntegralInfo> &integralInfos
                      ) {
  for (const auto &integral : integralInfos) {
    if ( ! a.isArgumentGiven(integral.name) ) continue;
    const auto& i(integral.indices);

    LOGGER(1) << "Computing " <<  integral.name << std::endl;

    std::vector<double> result;
    result = std::move(engine.compute(i[0], i[1], i[2], i[3]));
    std::vector<int> lens(4);
    for (unsigned int j(0) ; j < 4 ; j++) {
      lens[j] = (int)engine.indexToInt(i[j]);
    }
    a.allocatedTensorArgument<double>( integral.name
                                     , stdVectorToTensor(result, lens)
                                     );
  }
}

template<>
void computeAndExport ( Algorithm &a
                      , IntegralProvider< CTF::Tensor<double> > &engine
                      , std::vector<IntegralInfo> &integralInfos
                      ) {

  for (const auto &integral : integralInfos) {
    if ( ! a.isArgumentGiven(integral.name) ) continue;
    const auto& i(integral.indices);

    LOGGER(1) << "Computing " <<  integral.name << std::endl;

    a.allocatedTensorArgument<double>
      ( integral.name
      , new CTF::Tensor<double>(engine.compute(i[0], i[1], i[2], i[3]))
      );

  }
}

void CoulombIntegralsFromRotatedCoulombIntegrals::run() {

  std::vector<std::string> args( { "OrbitalCoefficients"
                                 , "CoulombIntegrals"
                                 , "nelec"
                                 , "No"
                                 , "chemistNotation"
                                 , "engine"
                                 } );

  auto C(getTensorArgument("OrbitalCoefficients"));
  auto V(getTensorArgument("CoulombIntegrals"));
  const int nelect(getIntegerArgument("nelec", -1));
  const int No(getIntegerArgument("No", nelect/2));
  const int Nv(V->lens[0] - No);
  const bool chemistNotation(getIntegerArgument("chemistNotation", 1) == 1);

  std::vector<IntegralInfo> integralInfos =
    { {"HHHHCoulombIntegrals", {NO,NO,NO,NO}, "ijkl"}
    , {"HHHPCoulombIntegrals", {NO,NO,NO,NV}, "ijka"}
    , {"HHPHCoulombIntegrals", {NO,NO,NV,NO}, "ijak"}
    , {"HHPPCoulombIntegrals", {NO,NO,NV,NV}, "ijab"}
    , {"HPHHCoulombIntegrals", {NO,NV,NO,NO}, "iajk"}
    , {"HPHPCoulombIntegrals", {NO,NV,NO,NV}, "iajb"}
    , {"HPPHCoulombIntegrals", {NO,NV,NV,NO}, "iabj"}
    , {"HPPPCoulombIntegrals", {NO,NV,NV,NV}, "iabc"}
    , {"PHHHCoulombIntegrals", {NV,NO,NO,NO}, "aijk"}
    , {"PHHPCoulombIntegrals", {NV,NO,NO,NV}, "aijb"}
    , {"PHPHCoulombIntegrals", {NV,NO,NV,NO}, "aibj"}
    , {"PHPPCoulombIntegrals", {NV,NO,NV,NV}, "aibc"}
    , {"PPHHCoulombIntegrals", {NV,NV,NO,NO}, "abij"}
    , {"PPHPCoulombIntegrals", {NV,NV,NO,NV}, "abic"}
    , {"PPPHCoulombIntegrals", {NV,NV,NV,NO}, "abci"}
    , {"PPPPCoulombIntegrals", {NV,NV,NV,NV}, "abcd"}
    };

  for (const auto& i: integralInfos) args.push_back(i.name);
  checkArgumentsOrDie(args);

  LOGGER(1) << "Note: CoulombIntegrals have to be in chemist notation!"
            << std::endl;
  LOGGER(1) << "Note: Phyiscs notation only implemented for ctf provider"
            << std::endl;
  LOGGER(1) << "No: " << No << std::endl;
  LOGGER(1) << "Nv: " << Nv << std::endl;


  struct EngineInfo {
    enum Name { CTF, VECTOR_SLOW, VECTOR_FAST };
    static Name fromString(const std::string name) {
      if (name == "ctf") return CTF;
      if (name == "VectorSlow") return VECTOR_SLOW;
      if (name == "VectorFast") return VECTOR_FAST;
      throw new EXCEPTION("Engine name not recognised");
    }
  };

  const EngineInfo::Name engineType(EngineInfo::fromString(
                                    getTextArgument("engine", "ctf")));

  LOGGER(1) << "Using engine: " << engineType << std::endl;
  LOGGER(1) << "Using chemistNotation?: " << chemistNotation << std::endl;

  if (engineType == EngineInfo::VECTOR_FAST) {
    VectorIntegralProvider engine(No, Nv, chemistNotation, *C, *V);
    computeAndExport(*this, engine, integralInfos);
  } else if (engineType == EngineInfo::VECTOR_SLOW) {
    SlowVectorIntegralProvider engine(No, Nv, chemistNotation, *C, *V);
    computeAndExport(*this, engine, integralInfos);
  } else {
    CtfIntegralProvider engine(No, Nv, chemistNotation, *C, *V);
    computeAndExport(*this, engine, integralInfos);
  }

  EMIT() << YAML::Key << "No" << YAML::Value << No
         << YAML::Key << "Nv" << YAML::Value << Nv
         << YAML::Key << "engine" << YAML::Value << engineType
         << YAML::Key << "chemist-notation" << YAML::Value << chemistNotation
         ;

}
