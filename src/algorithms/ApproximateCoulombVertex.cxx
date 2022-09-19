#include <algorithms/ApproximateCoulombVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(ApproximateCoulombVertex);

ApproximateCoulombVertex::ApproximateCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ApproximateCoulombVertex::~ApproximateCoulombVertex() {
}

void ApproximateCoulombVertex::run() {
  Tensor<complex> *GammaGqr(
    getTensorArgument<complex>("FullCoulombVertex")
  );
  GammaGqr->set_name("GammaGqr");
  Tensor<complex> *UGF(
    getTensorArgument<complex>("CoulombVertexSingularVectors")
  );
  Tensor<complex> UTGF(*UGF);
  conjugate(UTGF);
  int lens[] = { static_cast<int>(UGF->lens[1]),
                 static_cast<int>(GammaGqr->lens[1]),
                 static_cast<int>(GammaGqr->lens[2]) };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *GammaFqr = new Tensor<complex>(
    3, lens, syms, *GammaGqr->wrld, "GammaFqr"
  );
  allocatedTensorArgument<complex>(
    "CoulombVertex", GammaFqr
  );
  LOG(1, "ApproximateCoulombVertex")
    << "Using NF=" << UGF->lens[1] << " field variables to build " << GammaFqr->get_name()
    << " from " << GammaGqr->get_name() << " with NG=" << UGF->lens[0] << " grid points" << std::endl;
  (*GammaFqr)["Fqr"] = (*GammaGqr)["Gqr"] * UTGF["GF"];
}

void ApproximateCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex, DryTensor<complex>>("FullCoulombVertex")
  );
  DryTensor<complex> *UGF(
    getTensorArgument<complex, DryTensor<complex>>(
      "CoulombVertexSingularVectors"
    )
  );
  DryTensor<complex> UTGF(*UGF);
  int lens[] = { static_cast<int>(UGF->lens[1]),
                 static_cast<int>(GammaGqr->lens[1]),
                 static_cast<int>(GammaGqr->lens[2]) };
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *GammaFqr = new DryTensor<complex>(
    3, lens, syms, SOURCE_LOCATION
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "CoulombVertex", GammaFqr
  );
}

