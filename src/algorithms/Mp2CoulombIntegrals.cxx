#include <algorithms/Mp2CoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(Mp2CoulombIntegrals);

Mp2CoulombIntegrals::Mp2CoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2CoulombIntegrals::~Mp2CoulombIntegrals() {
}

/**
 * \brief Calculates Coulomb integrals Vabij from GammaGai Coulomb Vertex
 */
void Mp2CoulombIntegrals::run() {
  Tensor<complex> *GammaGpq(
    getTensorArgument<complex>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(
    getTensorArgument<>("HoleEigenEnergies")
  );
  Tensor<> *epsa(
    getTensorArgument<>("ParticleEigenEnergies")
  );

  LOG(0, "MP2CoulombIntegrals") <<
    "Reading Coulomb integrals form vertex " << GammaGpq->get_name() << " ...";

  // Compute the no,nv,nG,np
  int nG(GammaGpq->lens[0]);
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);
  int np = no + nv;

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { nv, nv, no, no };
  Tensor<> *Vabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij"));
  allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
  int GaiStart[] = {0 ,no, 0};
  int GaiEnd[]   = {nG,np,no};
  Tensor<complex> GammaGai(GammaGpq->slice(GaiStart,GaiEnd));

  // Split GammaGab,GammaGai,GammaGij into real and imaginary parts
  Tensor<> realGammaGai(
  3, GammaGai.lens, GammaGai.sym, *GammaGai.wrld, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
  3, GammaGai.lens, GammaGai.sym, *GammaGai.wrld, "ImagGammaGai"
  );
  fromComplexTensor(GammaGai, realGammaGai, imagGammaGai);

  // Compute the integrals Vabij
  (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
  (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];

  // Print okay
  LOG(0, "MP2CoulombIntegrals") << " OK" << std::endl;

  // Print test norm2 of GammaGai, GammaGab, GammaGij
  // GammaGai
  //double error(realGammaGai.norm2());
  //LOG(4) << "|realGammaGai| = " << error << std::endl;
  //error = imagGammaGai.norm2();
  //LOG(4) << "|imagGammaGai| = " << error << std::endl;
  // GammaGab
  //error = realGammaGab.norm2();
  //LOG(4) << "|realGammaGab| = " << error << std::endl;
  //error = imagGammaGab.norm2();
  //LOG(4) << "|imagGammaGab| = " << error << std::endl;
  // GammaGij
  //error = realGammaGij.norm2();
  //LOG(4) << "|realGammaGij| = " << error << std::endl;
  //error = imagGammaGij.norm2();
  //LOG(4) << "|imagGammaGij| = " << error << std::endl;

  // Print test norm2 of Vabij Vaibj Vaijb Vijkl Vabcd
  //error = Vabcd->norm2();
  //LOG(4) << "|Vabcd| = " << error << std::endl;
  //error = Vabij->norm2();
  //LOG(4) << "|Vabij| = " << error << std::endl;
  //error = Vaibj->norm2();
  //LOG(4) << "|Vaibj| = " << error << std::endl;
  //error = Vijkl->norm2();
  //LOG(4) << "|Vijkl| = " << error << std::endl;
}
