/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "FtodRankDecomposition.hpp"
#include "Exception.hpp"
#include "util/Log.hpp"
#include "util/ComplexPolynomialRootFinder.hpp"
#include "util/CubicPolynomialRootFinder.hpp"
#include <iostream>
#include <random>
#include <limits>


using namespace CTF;

FtodRankDecomposition::FtodRankDecomposition(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
}

FtodRankDecomposition::~FtodRankDecomposition() {
}

void FtodRankDecomposition::run() {
  rank = getIntegerArgument("rank");
  chiR = const_cast<Tensor<> *>(getTensorArgument("chiR"));
  chiI = const_cast<Tensor<> *>(getTensorArgument("chiI"));
  LOG(3) << "rank=" << rank << std::endl;
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  LOG(3) << "nG=" << nG << ", np=" << np << std::endl;
  // allocate decompisiton tensors
  X = new Matrix<>(int(rank), np, NS, *chiR->wrld, "XRp", chiR->profile);
  gamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamRRG", chiR->profile);
  gamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamIRG", chiR->profile);

  // allocate intermediate tensors:
  // tensor rank approximation
  chi0R = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0RGqr", chiR->profile);
  chi0I = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0IGqr", chiR->profile);
  // residuum
  RR = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RRGqr", chiR->profile);
  RI = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RIGqr", chiR->profile);

  // gradients
  dX = new Matrix<>(int(rank), np, NS, *chiR->wrld, "dXRp", chiR->profile);
  dGamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "dGamRRG", chiR->profile);
  dGamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "dGamIRG", chiR->profile);

  // NOTE that the elements don't need to be copied
  sX = new Matrix<>(*dX);
  sGamR = new Matrix<>(*dGamR);
  sGamI = new Matrix<>(*dGamI);

  initializeX();
  initializeGam();
  double epsilon = 1e-8;
  optimize(epsilon);
/*
  for (;;) {
    optimizeGam(epsilon);
    optimizeX(epsilon);
    double s(X->norm2());
    (*X)["Rp"] *= 1.0/s;
    (*gamR)["RG"] *= std::sqrt(s);
    (*gamI)["RG"] *= std::sqrt(s);
    epsilon /= 1000;
  }
*/
}


void FtodRankDecomposition::calculateChi0() {
  // TODO: try to keep this tensor allocated
  int np(chiR->lens[1]);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);

  x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
  (*chi0R)["Gqr"] = x["Rqr"]*(*gamR)["RG"];
  x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
  (*chi0I)["Gqr"] = x["Rqr"]*(*gamI)["RG"];
}

void FtodRankDecomposition::calculateResiduum() {
  // calculate the residuum
  (*RR)["Gqr"] = (*chi0R)["Gqr"] - (*chiR)["Gqr"];
  (*RI)["Gqr"] = (*chi0I)["Gqr"] - (*chiI)["Gqr"];
  Scalar<> s(*chiR->wrld);
  s[""] =  (*RR)["Gqr"] * (*RR)["Gqr"];
  s[""] += (*RI)["Gqr"] * (*RI)["Gqr"];
  R = s.get_val();
}

void FtodRankDecomposition::calculateGradient() {
  // TODO: try to keep these tensors allocated
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);
  int gLens[] = {int(rank), nG, np};
  int gSyms[] = {NS, NS, NS};
  Tensor<> g(3, gLens, gSyms, *chiR->wrld, "gRGq", chiR->profile);

  // calculate the gradient of X
  x["Rqr"] = (*RR)["Gqr"]*(*gamR)["RG"];
  // calculate 2*(x["Rpr"]*(*X)["Rr"] + x["Rpq"]*(*X)["Rq"])
  (*dX)["Rp"] = 4.0*x["Rpr"]*(*X)["Rr"];
  x["Rqr"] = (*RI)["Gqr"]*(*gamI)["RG"];
  (*dX)["Rp"] += 4.0*x["Rpr"]*(*X)["Rr"];

  // calculate the gradient of gamR and gamI
  g["RGq"] = (*RR)["Gqr"]*(*X)["Rr"];
  (*dGamR)["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
  g["RGq"] = (*RI)["Gqr"]*(*X)["Rr"];
  (*dGamI)["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
}

void FtodRankDecomposition::initializeRandom(Tensor<> &t, int64_t seed) {
  int64_t indicesCount, *indices;
  double *values;
  std::mt19937 random;
  random.seed(seed);
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  t.read_local(&indicesCount, &indices, &values);
  for (int64_t i(0); i < indicesCount; ++i) {
    values[i] = normalDistribution(random);
  }
  t.write(indicesCount, indices, values);
  free(indices); free(values);
}

void FtodRankDecomposition::initializeX() {
  initializeRandom(*X, chiR->wrld->rank);
}

void FtodRankDecomposition::initializeGam() {
  initializeRandom(*gamR, 2*chiR->wrld->rank+0);
  initializeRandom(*gamI, 2*chiR->wrld->rank+1);
}

void FtodRankDecomposition::lineSearchPart(
  Tensor<> &R, Tensor<> &G, Tensor<> &dG
) {
  int np(R.lens[1]);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  // TODO: try keeping these tensors for the
  // evaluation of the real and the imaginary part
  Tensor<> XX(3, xLens, xSyms, *R.wrld, "XXRqr", R.profile);
  Tensor<> dXX(3, xLens, xSyms, *R.wrld, "dXXRqr", R.profile);
  Tensor<> dXdX(3, xLens, xSyms, *R.wrld, "dXdXRqr", R.profile);
  XX["Rqr"] = (*X)["Rq"] * (*X)["Rr"];
  dXX["Rqr"]  =  (*sX)["Rq"] * (*X)["Rr"];
  dXX["Rqr"] +=  (*X)["Rq"] * (*sX)["Rr"];
  dXdX["Rqr"] = (*sX)["Rq"] * (*sX)["Rr"];

  // TODO: no copy needed for the following allocations
  // TODO: try keeping these tensors allocated for the
  // evaluation of the real and the imaginary part
  Tensor<> dXXG(R);
  Tensor<> XXdG(R);
  Tensor<> dXdXG(R);
  Tensor<> dXXdG(R);
  Tensor<> dXdXdG(R);
  // 1 alpha:
  dXXG["Gqr"] = dXX["Rqr"] * G["RG"];
  XXdG["Gqr"] = XX["Rqr"] * dG["RG"];
  // 2 alphas:
  dXdXG["Gqr"] = dXdX["Rqr"] * G["RG"];
  dXXdG["Gqr"] = dXX["Rqr"] * dG["RG"];
  // 3 alphas:
  dXdXdG["Gqr"] = dXdX["Rqr"] * dG["RG"];

  // alpha^0:
  Scalar<> s(*R.wrld);
  s[""] = R["Gqr"] * R["Gqr"];
  a[0] += s.get_val();

  // alpha^1:
  s[""] =  2.0 * dXXG["Gqr"] * R["Gqr"];
  s[""] += 2.0 * XXdG["Gqr"] * R["Gqr"];
  a[1] += s.get_val();

  // alpha^2:
  s[""] =  2.0 * dXdXG["Gqr"] * R["Gqr"];
  s[""] += 2.0 * dXXdG["Gqr"] * R["Gqr"];
  s[""] +=       dXXG["Gqr"] * dXXG["Gqr"];
  s[""] += 2.0 * dXXG["Gqr"] * XXdG["Gqr"];
  s[""] +=       XXdG["Gqr"] * XXdG["Gqr"];
  a[2] += s.get_val();

  // alpha^3:
  s[""] =  2.0 * dXdXdG["Gqr"] * R["Gqr"];
  s[""] += 2.0 * dXdXG["Gqr"] * dXXG["Gqr"];
  s[""] += 2.0 * dXdXG["Gqr"] * XXdG["Gqr"];
  s[""] += 2.0 * dXXdG["Gqr"] * dXXG["Gqr"];
  s[""] += 2.0 * dXXdG["Gqr"] * XXdG["Gqr"];
  a[3] += s.get_val();

  // alpha^4:
  s[""] =  2.0 * dXdXdG["Gqr"] * dXXG["Gqr"];
  s[""] += 2.0 * dXdXdG["Gqr"] * XXdG["Gqr"];
  s[""] +=       dXdXG["Gqr"] * dXdXG["Gqr"];
  s[""] += 2.0 * dXdXG["Gqr"] * dXXdG["Gqr"];
  s[""] +=       dXXdG["Gqr"] * dXXdG["Gqr"];
  a[4] += s.get_val();

  // alpha^5:
  s[""] =  2.0 * dXdXdG["Gqr"] * dXdXG["Gqr"];
  s[""] += 2.0 * dXdXdG["Gqr"] * dXXdG["Gqr"];
  a[5] += s.get_val();
  // alpha^6:
  s[""] = dXdXdG["Gqr"] * dXdXdG["Gqr"];
  a[6] += s.get_val();
}

double FtodRankDecomposition::lineSearch() {
  for (int k(0); k <= 6; ++k) {
    a[k] = 0.0;
  }
  lineSearchPart(*RR, *gamR, *sGamR);
  lineSearchPart(*RI, *gamI, *sGamI);
  util::Polynomial<double> p(a, 6);
  // the derivative polynomial
  // TODO: add retrieving derivative polynomial to Polynomial class
  // TODO: add element cast to Polynomial class
  for (int k(0); k <= 6; ++k) LOG(4) << "a[" << k << "]=" << a[k] << std::endl;
  std::complex<double> b[6];
  for (int k(1); k <= 6; ++k) {
    b[k-1] = k*a[k];
  }
  for (int k(0); k <= 5; ++k) LOG(4) << "b[" << k << "]=" << b[k] << std::endl;
  util::ScaledPolynomial<std::complex<double>> dp(b, 5);
  util::ComplexPolynomialRootFinder rootFinder(&dp);
  rootFinder.findRoots();
  double min(std::numeric_limits<double>::infinity());
  double argmin(0.0);
  for (int i(0); i < 6; ++i) {
    std::complex<double> root = rootFinder.getRoot(i);
    // FIXME: control accuarcy:
    if (std::abs(root.imag()) < 1e-10) {
      root /= dp.getScale();
      double r(p.at(root.real()));
      if (r < min) {
        argmin = root.real();
        min = r;
      }
    }
  }
  
  LOG(3) << "alpha=" << argmin << std::endl;
  LOG(3) << "min(R)=" << min << std::endl;
  return argmin;
}

void FtodRankDecomposition::optimize(double const epsilon) {
  calculateChi0();
  calculateResiduum();
  calculateGradient();
  (*sX)["Rp"] = -1.0*(*dX)["Rp"];
  (*sGamR)["RG"] = -1.0*(*dGamR)["RG"];
  (*sGamI)["RG"] = -1.0*(*dGamI)["RG"];
  double alpha(lineSearch());
  (*X)["Rp"] += alpha*(*sX)["Rp"];
  (*gamR)["RG"] += alpha*(*sGamR)["RG"];
  (*gamI)["RG"] += alpha*(*sGamI)["RG"];
  for (int count(0); ; ++count) {
    Tensor<> lastDX(*dX);
    Tensor<> lastDGamR(*dGamR);
    Tensor<> lastDGamI(*dGamI);
    calculateChi0();
    calculateResiduum();
    LOG(2) << count << ":" << std::endl;
    LOG(2) << "  R=" << R << std::endl;
    calculateGradient();
    Scalar<> s(*lastDX.wrld);
    s[""] = lastDX["Rp"] * lastDX["Rp"];
    s[""] += lastDGamR["RG"] * lastDGamR["RG"];
    s[""] += lastDGamI["RG"] * lastDGamI["RG"];
    double beta(s.get_val());
    LOG(2) << "  |lastD|=" << beta << std::endl;
    lastDX["Rp"] -= (*dX)["Rp"];
    lastDGamR["RG"] -= (*dGamR)["RG"];
    lastDGamI["RG"] -= (*dGamI)["RG"];
    s[""] = (*dX)["Rp"] * lastDX["Rp"];
    s[""] += (*dGamR)["Rp"] * lastDGamR["Rp"];
    s[""] += (*dGamI)["Rp"] * lastDGamI["Rp"];
    beta = std::max(0.0, -s.get_val() / beta);
    LOG(2) << "  beta=" << beta << std::endl;
    (*sX)["Rp"] = beta*(*sX)["Rp"] - (*dX)["Rp"];
    (*sGamR)["RG"] = beta*(*sGamR)["RG"] - (*dGamR)["RG"];
    (*sGamI)["RG"] = beta*(*sGamI)["RG"] - (*dGamI)["RG"];
    double alpha(lineSearch());
    if (alpha < epsilon) return;
    (*X)["Rp"] += alpha*(*sX)["Rp"];
    (*gamR)["RG"] += alpha*(*sGamR)["RG"];
    (*gamI)["RG"] += alpha*(*sGamI)["RG"];
  }
}


void FtodRankDecomposition::lineSearchXPart(
  Tensor<> &chi0, Tensor<> &chi, Tensor<> &gam
) {
  // TODO: try to keep these tensors allocated
  int np(chi.lens[1]);
  // TODO: no copy needed for the following two allocations
  Tensor<> dChi0(chi);
  Tensor<> ddChi0(chi);
  Tensor<> deltaChi0(chi0);
  deltaChi0["Gqr"] -= chi["Gqr"];
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chi.wrld, "xRqr", chi.profile);
  x["Rqr"] = (*sX)["Rq"] * (*X)["Rr"];
  dChi0["Gqr"] = x["Rqr"] * gam["RG"];
//  dChi0["Gqr"] += x["Rrq"] * gam["RG"];
  // FIXME: check if tensors can be used in update
  dChi0["Gqr"] += dChi0["Grq"];
  x["Rqr"] = (*sX)["Rq"] * (*sX)["Rr"];
  ddChi0["Gqr"] = x["Rqr"] * gam["RG"];
  Scalar<> s(*chi.wrld);
  s[""] = deltaChi0["Gqr"] * deltaChi0["Gqr"];
  a0 += s.get_val();
  s[""] = 2.0*dChi0["Gqr"] * deltaChi0["Gqr"];
  a1 += s.get_val();
  s[""] = 2.0*ddChi0["Gqr"] * deltaChi0["Gqr"];
  s[""] += dChi0["Gqr"] * dChi0["Gqr"];
  a2 += s.get_val();
  s[""] = 2.0*ddChi0["Gqr"] * dChi0["Gqr"];
  a3 += s.get_val();
  s[""] = ddChi0["Gqr"] * ddChi0["Gqr"];
  a4 += s.get_val();
}

double FtodRankDecomposition::lineSearchX() {
  a0 = a1 = a2 = a3 = a4 = 0.0;
  lineSearchXPart(*chi0R, *chiR, *gamR);
  lineSearchXPart(*chi0I, *chiI, *gamI);
  util::CubicPolynomialRootFinder rootFinder(
    a1, 2.0*a2, 3.0*a3, 4.0*a4
  );
  double roots[3];
  int rootsCount;
  rootsCount = rootFinder.findRoots(roots);
  double min(std::numeric_limits<double>::infinity());
  double argmin(0.0);
  for (int i(0); i < rootsCount; ++i) {
    if (rootFinder.evaluateDerivativeAt(roots[i]) > 0.0) {
      // only if second derivative is positive consider it
      double r(a4*roots[i]+a3);
      r = r*roots[i] + a2;
      r = r*roots[i] + a1;
      r = r*roots[i] + a0;
      if (r < min) {
        argmin = roots[i];
        min = r;
      }
    }
  }
  LOG(3) << "alpha=" << argmin << std::endl;
  LOG(3) << "min(R)=" << min << std::endl;
  return argmin;
}

void FtodRankDecomposition::optimizeX(double const epsilon) {
  calculateChi0();
  calculateResiduum();
  calculateGradient();
  (*sX)["Rp"] = -1.0*(*dX)["Rp"];
  double alpha(lineSearchX());
  (*X)["Rp"] += alpha*(*sX)["Rp"];
  for (int count(0); ; ++count) {
    Tensor<> lastDX(*dX);
    calculateChi0();
    calculateResiduum();
    LOG(2) << "X " << count << ":" << std::endl;
    LOG(2) << "  R=" << R << std::endl;
    calculateGradient();
    Scalar<> s(*lastDX.wrld);
    s[""] = lastDX["Rp"] * lastDX["Rp"];
    double beta(s.get_val());
    LOG(2) << "  |lastDX|=" << beta << std::endl;
    lastDX["Rp"] -= (*dX)["Rp"];
    s[""] = (*dX)["Rp"] * lastDX["Rp"];
    beta = std::max(0.0, -s.get_val() / beta);
    LOG(2) << "  beta=" << beta << std::endl;
    (*sX)["Rp"] = beta*(*sX)["Rp"] - (*dX)["Rp"];
    double alpha(lineSearchX());
    if (alpha < epsilon) return;
    (*X)["Rp"] += alpha*(*sX)["Rp"];
  }
}

void FtodRankDecomposition::lineSearchGamPart(
  Tensor<> &chi0, Tensor<> &chi, Tensor<> &sGam
) {
  // TODO: try to keep these tensors allocated
  int np(chi.lens[1]);
  // TODO: no copy needed for the following allocation
  Tensor<> dChi0(chi);
  Tensor<> deltaChi0(chi0);
  deltaChi0["Gqr"] -= chi["Gqr"];
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chi.wrld, "xRqr", chi.profile);
  x["Rqr"] = (*X)["Rq"] * (*X)["Rr"];
  dChi0["Gqr"] = x["Rqr"] * sGam["RG"];
  Scalar<> s(*chi.wrld);
  s[""] = deltaChi0["Gqr"] * deltaChi0["Gqr"];
  a0 += s.get_val();
  s[""] = 2.0*dChi0["Gqr"] * deltaChi0["Gqr"];
  a1 += s.get_val();
  s[""] = dChi0["Gqr"] * dChi0["Gqr"];
  a2 += s.get_val();
}


double FtodRankDecomposition::lineSearchGam() {
  a0 = a1 = a2 = a3 = a4 = 0.0;
  lineSearchGamPart(*chi0R, *chiR, *sGamR);
  lineSearchGamPart(*chi0I, *chiI, *sGamI);
  double argmin(-0.5*a1/a2);
  double min(a2*argmin + a1);
  min = min*argmin + a0;
  LOG(3) << "alpha=" << argmin << std::endl;
  LOG(3) << "min(R)=" << min << std::endl;
  return argmin;
}

void FtodRankDecomposition::optimizeGam(double const epsilon) {
  calculateChi0();
  calculateResiduum();
  calculateGradient();
  (*sGamR)["RG"] = -1.0*(*dGamR)["RG"];
  (*sGamI)["RG"] = -1.0*(*dGamI)["RG"];
  double alpha(lineSearchGam());
  (*gamR)["RG"] += alpha*(*sGamR)["RG"];
  (*gamI)["RG"] += alpha*(*sGamI)["RG"];
  for (int count(0); ; ++count) {
    Tensor<> lastDGamR(*dGamR);
    Tensor<> lastDGamI(*dGamI);
    calculateChi0();
    calculateResiduum();
    LOG(2) << "G " << count << ":" << std::endl;
    LOG(2) << "  R=" << R << std::endl;
    calculateGradient();
    Scalar<> s(*lastDGamR.wrld);
    s[""] = lastDGamR["RG"] * lastDGamR["RG"];
    s[""] += lastDGamI["RG"] * lastDGamI["RG"];
    double beta(s.get_val());
    LOG(2) << "  |lastDGam|=" << beta << std::endl;
    lastDGamR["RG"] -= (*dGamR)["RG"];
    lastDGamI["RG"] -= (*dGamI)["RG"];
    s[""] = (*dGamR)["RG"] * lastDGamR["RG"];
    s[""] += (*dGamI)["RG"] * lastDGamI["RG"];
    beta = std::max(0.0, -s.get_val() / beta);
    LOG(2) << "  beta=" << beta << std::endl;
    (*sGamR)["RG"] = beta*(*sGamR)["RG"] - (*dGamR)["RG"];
    (*sGamI)["RG"] = beta*(*sGamI)["RG"] - (*dGamI)["RG"];
    double alpha(lineSearchGam());
    if (alpha < epsilon) return;
    (*gamR)["RG"] += alpha*(*sGamR)["RG"];
    (*gamI)["RG"] += alpha*(*sGamI)["RG"];
  }
}

void FtodRankDecomposition::testGradient() {
  calculateChi0();
  calculateResiduum();
  calculateGradient();
  double R0(R);
  Matrix<> oldX(*X), oldGamR(*gamR), oldGamI(*gamI);
  double delta(1e-10);
  double accuracy(1e-6);
//  for (double delta(width/100); delta <= width; delta += width/100) {
    (*X)["Rq"] = oldX["Rq"]+delta*(*dX)["Rq"];
    (*gamR)["RG"] = oldGamR["RG"]+delta*(*dGamR)["RG"];
    (*gamI)["RG"] = oldGamI["RG"]+delta*(*dGamI)["RG"];
    calculateChi0();
    calculateResiduum();
    double epsilon(R - R0);
    Scalar<> s(*chiR->wrld);
    s[""] =  (*dX)["Rq"] * (*dX)["Rq"];
    s[""] += (*dGamR)["RG"] * (*dGamR)["RG"];
    s[""] += (*dGamI)["RG"] * (*dGamI)["RG"];
    double dotEpsilon(delta * s.get_val());
    if (std::abs(dotEpsilon-epsilon)/epsilon > accuracy)
      throw new Exception("wrong gradient");
//  }
}
