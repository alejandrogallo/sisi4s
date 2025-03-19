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
#define HANDLE_INTEGRAL(indices, no, nv)                                       \
  do {                                                                         \
    if (in.present(#indices)) {                                                \
      auto _v = in.get<Tensor<double> *>(#indices);                            \
      V->with_##indices(_v);                                                   \
      if (no > 0) { V->No = _v->lens[no]; }                                    \
      if (nv > 0) { V->Nv = _v->lens[nv]; }                                    \
    }                                                                          \
  } while (0)

  HANDLE_INTEGRAL(hhhh, 0, -1);
  HANDLE_INTEGRAL(hhhp, 0, 3);
  HANDLE_INTEGRAL(hhph, 0, 2);
  HANDLE_INTEGRAL(hhpp, 0, 2);
  HANDLE_INTEGRAL(hphh, 0, 1);
  HANDLE_INTEGRAL(hphp, 0, 1);
  HANDLE_INTEGRAL(hpph, 0, 1);
  HANDLE_INTEGRAL(hppp, 0, 1);
  HANDLE_INTEGRAL(phhh, 1, 0);
  HANDLE_INTEGRAL(phhp, 1, 0);
  HANDLE_INTEGRAL(phph, 1, 0);
  HANDLE_INTEGRAL(phpp, 1, 0);
  HANDLE_INTEGRAL(pphh, 2, 0);
  HANDLE_INTEGRAL(pphp, 2, 0);
  HANDLE_INTEGRAL(ppph, 3, 0);
  HANDLE_INTEGRAL(pppp, -1, 1);

  out.set<sisi4s::CoulombIntegrals<double> *>("out", V);
  LOG(0, "CoulombIntegrals") << "No: " << V->No << std::endl;
  LOG(0, "CoulombIntegrals") << "Nv: " << V->Nv << std::endl;
}

} // namespace step
} // namespace sisi4s
