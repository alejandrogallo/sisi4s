#include <algorithms/DoublesAmplitudesDecomposition.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <DryTensor.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;
using std::make_shared;

ALGORITHM_REGISTRAR_DEFINITION(DoublesAmplitudesDecomposition);

DoublesAmplitudesDecomposition::DoublesAmplitudesDecomposition(
    std::vector<Argument> const &argumentList)
    : Algorithm(argumentList) {}

DoublesAmplitudesDecomposition::~DoublesAmplitudesDecomposition() {}

void DoublesAmplitudesDecomposition::run() {
  diagonlizeAmplitudes();
  sliceLargestEigenValues();
}

void DoublesAmplitudesDecomposition::dryRun() {
  // Read the Coulomb vertex GammaGai
  DryTensor<> *Tabij(
      getTensorArgument<double, DryTensor<>>("DoublesAmplitudes"));

  Nv = Tabij->lens[0] * Tabij->lens[1];
  // get number of field variables
  int NL(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NL == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NL = static_cast<int>(Nv * reduction + 0.5);
  }
}

void DoublesAmplitudesDecomposition::diagonlizeAmplitudes() {
  Tensor<double> *Tabij(getTensorArgument<>("DoublesAmplitudes"));

  // reorder to Taibj
  Nv = Tabij->lens[0];
  No = Tabij->lens[2];
  int lens[] = {Nv, No, Nv, No};
  Taibj =
      make_shared<Tensor<double>>(4, lens, Tabij->sym, *Tabij->wrld, "Taibj");
  (*Taibj)["aibj"] = (*Tabij)["abij"];

  // T(ai)(bj) = U.Lambda.U^T, seen as a matrix with compound indices
  BlacsWorld world(Taibj->wrld->rank, Taibj->wrld->np);
  NvNo = Nv * No;
  int scaTLens[2] = {NvNo, NvNo};
  auto scaTaibj(make_shared<ScaLapackMatrix<>>(*Taibj, scaTLens, &world));
  // release unneeded resources early
  Taibj.reset();

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(make_shared<ScaLapackMatrix<>>(*scaTaibj));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTaibj, scaU);
  lambdas = new double[NvNo];
  eigenSystem.solve(lambdas);
  scaTaibj.reset();

  // write matrix U(ai)(F) back to CTF as tensor UaiF
  int ULens[3] = {Nv, No, NvNo};
  UaiF =
      make_shared<Tensor<double>>(3, ULens, Tabij->sym, *Tabij->wrld, "UaiF");
  scaU->write(*UaiF);
  scaU.reset();

  // write Lambda and conj(sqrt(Lambda)) back to CTF
  lambdasCount = UaiF->wrld->rank == 0 ? NvNo : 0;
  lambdaIndices = new int64_t[lambdasCount];
  sqrtLambdas = new complex[lambdasCount];
  for (int64_t i(0); i < lambdasCount; ++i) {
    lambdaIndices[i] = i;
    sqrtLambdas[i] = sqrt(complex(lambdas[i]));
  }
  sqrtLambdaF = make_shared<Tensor<complex>>(1, &NvNo, UaiF->sym, *UaiF->wrld);
  sqrtLambdaF->set_name("sqrtLambda");
  sqrtLambdaF->write(lambdasCount, lambdaIndices, sqrtLambdas);
  LambdaF = new Tensor<double>(1, &NvNo, UaiF->sym, *UaiF->wrld);
  LambdaF->set_name("Lambda");
  LambdaF->write(lambdasCount, lambdaIndices, lambdas);
  delete[] lambdaIndices;
  delete[] sqrtLambdas;

  allocatedTensorArgument<>("DoublesAmplitudesEigenValues", LambdaF);
}

void DoublesAmplitudesDecomposition::sliceLargestEigenValues() {
  // get number of field variables
  NF = getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES);
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(Nv * reduction + 0.5);
  }
  NF = std::min(std::max(0, NF), NvNo);

  // identify largest NF |eigenvalues|
  int lower(0), upper(NvNo - 1);
  while (NvNo - 1 - upper + lower < NF) {
    if (std::abs(lambdas[lower]) > std::abs(lambdas[upper])) {
      ++lower;
    } else {
      --upper;
    }
  }
  delete[] lambdas;
  LOG(1, "DoublesAmplitudesDecomposition")
      << "taken values: 0<=i<" << lower << " and " << upper + 1 << "<=i<"
      << NvNo << std::endl;

  // slice eigen vectors U and values Lambda corresponding to NF largest values
  int LaiFLens[] = {Nv, No, NF};
  auto LaiF(
      make_shared<Tensor<double>>(3, LaiFLens, UaiF->sym, *UaiF->wrld, "LaiF"));
  auto LF(make_shared<Tensor<complex>>(1, &NF, UaiF->sym, *UaiF->wrld, "LF"));
  {
    int lowerStart[] = {0, 0, 0};
    int lowerEnd[] = {Nv, No, lower};
    LaiF->slice(lowerStart, lowerEnd, 0.0, *UaiF, lowerStart, lowerEnd, 1.0);
    LF->slice(&lowerStart[2],
              &lowerEnd[2],
              0.0,
              *sqrtLambdaF,
              &lowerStart[2],
              &lowerEnd[2],
              1.0);
    int upperStart[] = {0, 0, upper + 1};
    int upperEnd[] = {Nv, No, NvNo};
    int targetStart[] = {0, 0, lower};
    LaiF->slice(targetStart, LaiFLens, 0.0, *UaiF, upperStart, upperEnd, 1.0);
    LF->slice(&targetStart[2],
              &LaiFLens[2],
              0.0,
              *sqrtLambdaF,
              &upperStart[2],
              &upperEnd[2],
              1.0);
  }
  UaiF.reset();
  sqrtLambdaF.reset();

  auto cLaiF(make_shared<Tensor<complex>>(3,
                                          LaiFLens,
                                          LaiF->sym,
                                          *LaiF->wrld,
                                          "cLaiF"));
  toComplexTensor(*LaiF, *cLaiF);
  LaiF.reset();

  int LFaiLens[] = {NF, Nv, No};
  auto LFai(new Tensor<complex>(3, LFaiLens, cLaiF->sym, *cLaiF->wrld, "LFai"));
  (*LFai)["Fai"] = (*cLaiF)["aiF"] * (*LF)["F"];
  cLaiF.reset();
  LF.reset();

  allocatedTensorArgument<complex>("DoublesAmplitudesVertex", LFai);

  // recompose for testing
  Tensor<double> *Tabij(getTensorArgument<>("DoublesAmplitudes"));
  double tNorm(frobeniusNorm(*Tabij));
  Tensor<complex> cTabij(4, Tabij->lens, Tabij->sym, *Tabij->wrld, "cTabij");
  toComplexTensor(*Tabij, cTabij);
  cTabij["abij"] += (-1.0) * (*LFai)["Fai"] * (*LFai)["Fbj"];
  double dNorm(frobeniusNorm(cTabij));
  LOG(1, "DoublesAmplitudesDecomposition")
      << "|T-L.L|/|T|=" << dNorm / tNorm << std::endl;
}
