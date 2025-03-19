#include <Step.hpp>
#include <equations/CoulombIntegrals.hpp>

#define COULOMB_INTEGRAL(indices)                                              \
  {                                                                            \
#    indices, SPEC_VARIN(#indices " part of integrals", Tensor <double> *)     \
  }

namespace sisi4s {

namespace step {

DEFSPEC(CoulombIntegrals,
        SPEC_IN(COULOMB_INTEGRAL(hhhh),
                COULOMB_INTEGRAL(hhhp),
                COULOMB_INTEGRAL(hhph),
                COULOMB_INTEGRAL(hhpp),
                COULOMB_INTEGRAL(hphh),
                COULOMB_INTEGRAL(hphp),
                COULOMB_INTEGRAL(hpph),
                COULOMB_INTEGRAL(hppp),
                COULOMB_INTEGRAL(phhh),
                COULOMB_INTEGRAL(phhp),
                COULOMB_INTEGRAL(phph),
                COULOMB_INTEGRAL(phpp),
                COULOMB_INTEGRAL(pphh),
                COULOMB_INTEGRAL(pphp),
                COULOMB_INTEGRAL(ppph),
                COULOMB_INTEGRAL(pppp), ),
        SPEC_OUT({"out",
                  SPEC_VAROUT("TODO: DOC", sisi4s::CoulombIntegrals<double> *)
                      ->require()}, ));

DEFSTEP(CoulombIntegrals) {
  auto V = new sisi4s::CoulombIntegrals<double>();
#define HANDLE_INTEGRAL(indices)                                               \
  do {                                                                         \
    if (in.present(#indices)) {                                                \
      auto _v = in.get<Tensor<double> *>(#indices);                            \
      V->with_##indices(_v);                                                   \
    }                                                                          \
  } while (0)

  HANDLE_INTEGRAL(hhhh);
  HANDLE_INTEGRAL(hhhp);
  HANDLE_INTEGRAL(hhph);
  HANDLE_INTEGRAL(hhpp);
  HANDLE_INTEGRAL(hphh);
  HANDLE_INTEGRAL(hphp);
  HANDLE_INTEGRAL(hpph);
  HANDLE_INTEGRAL(hppp);
  HANDLE_INTEGRAL(phhh);
  HANDLE_INTEGRAL(phhp);
  HANDLE_INTEGRAL(phph);
  HANDLE_INTEGRAL(phpp);
  HANDLE_INTEGRAL(pphh);
  HANDLE_INTEGRAL(pphp);
  HANDLE_INTEGRAL(ppph);
  HANDLE_INTEGRAL(pppp);

  out.set<sisi4s::CoulombIntegrals<double> *>("out", V);
  LOG(0, "CoulombIntegrals") << "No: " << V->No << std::endl;
  LOG(0, "CoulombIntegrals") << "Nv: " << V->Nv << std::endl;
}

} // namespace step
} // namespace sisi4s
