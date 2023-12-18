#include <algorithms/BasisSetExtrapolation.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <DryTensor.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <util/Tensor.hpp>
#include <util/MpiCommunicator.hpp>
#include <math/PseudoInverseSvd.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(BasisSetExtrapolation) {}

DEFSPEC(
    BasisSetExtrapolation,
    SPEC_IN({"gamma", SPEC_VALUE_DEF("TODO: DOC", double, 1)},
            {"maxG", SPEC_VALUE_DEF("TODO: DOC", double, -1)},
            {"minG", SPEC_VALUE_DEF("TODO: DOC", double, -1)},
            {"volume", SPEC_POSITIVE("TODO: DOC", double)->require()},
            {"calculateQGG", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
            {"fitF12",
             SPEC_ONE_OF("It should be one or 2", int64_t, 1, 2)->require()},
            {"invertQGG", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
            {"iterations", SPEC_VALUE_DEF("TODO: DOC", int64_t, 10)},
            {"orbitalPairEnd", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},
            {"orbitalPairStart", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},
            {"slice", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},
            {"CoulombVertex",
             SPEC_VARIN("TODO: DOC", Tensor<complex> *)->require()},
            {"ParticleHoleCoulombVertex",
             SPEC_VARIN("TODO: DOC", Tensor<complex> *)},
            {"CoulombKernel", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"HoleEigenEnergies",
             SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},
            {"Momenta", SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},
            {"QGG", SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()},
            {"StructureFactor",
             SPEC_VARIN("TODO: DOC", Tensor<double> *)->require()}),
    SPEC_OUT({"f12", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"FittedSF", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"QGGd", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"QGGs", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"QGGt", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"QGGx", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(BasisSetExtrapolation) {

  bool fQGG(in.get<bool>("calculateQGG"));
  if (fQGG) {
    // slice QGG evalution in case of memory bottleneck for the exchange term
    int slice(in.get<int64_t>("slice"));
    int orbitalPairStart(in.get<int64_t>("orbitalPairStart"));
    int orbitalPairEnd(in.get<int64_t>("orbitalPairEnd"));
    LOG(0, "BasisSetExtrapolation:") << "evaluating QGG" << std::endl;
    evaluateQGG(orbitalPairStart, orbitalPairEnd, slice);
  }

  int fFitF12(in.get<int64_t>("fitF12"));
  if (fFitF12 == 1 || fFitF12 == 2) {
    real minG(in.get<double>("minG"));
    real maxG(in.get<double>("maxG"));
    if (minG > maxG || minG < 0.)
      throw new EXCEPTION("need fitting range:minG and maxG");
    LOG(0, "BasisSetExtrapolation") << "fitting gamma" << std::endl;
    fitF12(fFitF12, minG, maxG);
  }

  if (in.get<bool>("invertQGG")) {
    LOG(0, "BasisSetExtrapolation")
        << "inverting Q(G,G') --> correlation Factor" << std::endl;
    invertQGG();
  }
}

void BasisSetExtrapolation::evaluateQGG(int orbitalPairStart,
                                        int orbitalPairEnd,
                                        int slice) {

  PTR(Tensor<complex>) GammaGai;

  if (in.present("ParticleHoleCoulombVertex")) {
    GammaGai = NEW(Tensor<complex>,
                   in.get<Tensor<complex> *>("ParticleHoleCoulombVertex"));
  } else if (in.present("CoulombVertex")) {

    auto epsi =
        NEW(Tensor<double>, in.get<Tensor<double> *>("HoleEigenEnergies"));
    int No(epsi->lens[0]);
    auto GammaGqr =
        NEW(Tensor<complex>, in.get<Tensor<complex> *>("CoulombVertex"));
    int NG(GammaGqr->lens[0]);
    int Np(GammaGqr->lens[1]);
    int Nv(Np - No);
    int aStart(Np - Nv), aEnd(Np);
    int iStart(0), iEnd(No);
    int GaiStart[] = {0, aStart, iStart};
    int GaiEnd[] = {NG, aEnd, iEnd};
    GammaGai = NEW(Tensor<complex>, GammaGqr->slice(GaiStart, GaiEnd));
  } else {
    throw new EXCEPTION("Need Appropriate Coulomb Vertex");
  }

  // Gamma Gai --> Cai
  int NF(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);

  Tensor<double> *ctfCoulombKernel(in.get<Tensor<double> *>("CoulombKernel"));

  int NFF[] = {NF};
  Tensor<complex> invSqrtVG(1, NFF);

  CTF::Transform<real, complex>(
      std::function<void(real, complex &)>([](real vG, complex &invVG) {
        if (std::isinf(vG)) {
          invVG = 0.;
        } else {
          invVG = std::sqrt(1. / vG);
        }
      }))((*ctfCoulombKernel)["G"], invSqrtVG["G"]);

  (*GammaGai)["Gai"] *= invSqrtVG["G"];

  // determine if we are using full or half mesh

  auto momenta(NEW(Tensor<double>, in.get<Tensor<double> *>("Momenta")));
  sisi4s::Vector<> *cartesianMomenta(new sisi4s::Vector<>[NF]);
  momenta->read_all(&cartesianMomenta[0][0]);

  sisi4s::Vector<> check_grid;
  for (int g(0); g < NF; ++g) { check_grid += cartesianMomenta[g]; }
  PTR(Tensor<complex>) CGai;
  if (check_grid.length() > 1e-5) {
    // maybe the easiest is to double Cia(G).
    // we need to rescale Cia(G)=>Cia(G)/sqrt(2)
    LOG(1, "Build up Q(G,F)") << "working with half mesh" << std::endl;
    int NGai[] = {2 * NF - 1, Nv, No};
    real invsqrt(1. / std::sqrt(2.));
    CGai = NEW(Tensor<complex>, 3, NGai);
    // First put all 'positive' G in the full Cia(G)
    CGai->slice(std::vector<int>({0, 0, 0}).data(),
                std::vector<int>({NF, Nv, No}).data(),
                1.0,
                *GammaGai,
                std::vector<int>({0, 0, 0}).data(),
                std::vector<int>({NF, Nv, No}).data(),
                invsqrt);
    // Now all 'negative' G
    conjugate(*GammaGai);
    CGai->slice(std::vector<int>({NF, 0, 0}).data(),
                std::vector<int>({2 * NF - 1, Nv, No}).data(),
                1.0,
                *GammaGai,
                std::vector<int>({1, 0, 0}).data(),
                std::vector<int>({NF, Nv, No}).data(),
                invsqrt);
    NF = NF * 2 - 1;
  } else {
    LOG(1, "Build up Q(G,F)") << "working with full mesh" << std::endl;
    CGai = NEW(Tensor<complex>, *GammaGai);
  }

  if (orbitalPairStart < 0 || orbitalPairStart > No) { orbitalPairStart = 0; }
  if (orbitalPairEnd <= orbitalPairStart || orbitalPairEnd > No) {
    orbitalPairEnd = No;
  }

  LOG(1, "Orbital Pair analysis")
      << "Considering electron pairs " << orbitalPairStart << " to "
      << orbitalPairEnd << std::endl;
  No = orbitalPairEnd - orbitalPairStart;

  int sCGaiStart[] = {0, 0, orbitalPairStart};
  int sCGaiEnd[] = {NF, Nv, orbitalPairEnd};

  auto sCGai(NEW(Tensor<complex>, CGai->slice(sCGaiStart, sCGaiEnd)));

  auto conjCGai(NEW(Tensor<complex>, *sCGai));
  conjugate(*conjCGai);

  int NGG[] = {NF, NF};
  auto QGGs(new Tensor<complex>(2, NGG));
  auto QGGt(new Tensor<complex>(2, NGG));

  // Direct part.

  PTR(Tensor<complex>) FGone;
  PTR(Tensor<complex>) FGtwo;
  FGone = NEW(Tensor<complex>, 2, NGG);
  FGtwo = NEW(Tensor<complex>, 2, NGG);

  // oldschool
  //   (*FGone)["GF"] =  (*conjCGai)["Gai"] * (*sCGai)["Fai"];
  //   (*FGtwo)["GF"] =  (*conjCGai)["Fbj"] * (*sCGai)["Gbj"];
  //   (*QGGs)["GF"]  = (0.5) * (*FGone)["GF"] * (*FGtwo)["GF"];
  //   (*QGGt)["GF"]  = (1.5) * (*QGGs)["GF"];
  // newschool
  (*FGone)["GF"] = (*conjCGai)["Gbj"] * (*conjCGai)["Fbj"];
  (*FGtwo)["GF"] = (*sCGai)["Gai"] * (*sCGai)["Fai"];
  (*QGGs)["GF"] = (0.5) * (*FGone)["GF"] * (*FGtwo)["GF"];
  (*QGGt)["GF"] = (1.5) * (*QGGs)["GF"];

  if (out.present("QGGd")) {
    auto QGGd(new Tensor<double>(2, NGG));
    fromComplexTensor(*QGGs, *QGGd);
    out.set<Tensor<double> *>("QGGd", QGGd);
  }

  // Exchange part, slicing.

  if (slice < 0 || slice > No) slice = No;

  LOG(0, "Number of bands treated simultaniously:") << slice << std::endl;

  int numberSlices(
      std::ceil(static_cast<double>(No) / static_cast<double>(slice)));
  PTR(Tensor<complex>) cQGGx;
  if (out.present("QGGx")) { cQGGx = NEW(Tensor<complex>, 2, NGG); }
  for (int ii(0); ii < numberSlices; ++ii) {

    int startBandSlice(ii * slice);
    int endBandSlice((ii + 1) * slice);
    endBandSlice = std::min(endBandSlice, No);
    int NFG[] = {NF, NF, No, endBandSlice - startBandSlice};
    LOG(0, "Slicing for memory reasons:")
        << ii + 1 << " From: " << startBandSlice << " to: " << endBandSlice
        << std::endl;

    PTR(Tensor<complex>) FGij;
    PTR(Tensor<complex>) FGji;

    FGij = NEW(Tensor<complex>, 4, NFG);
    FGji = NEW(Tensor<complex>, 4, NFG);

    int CGajStart[] = {0, 0, startBandSlice};
    int CGajEnd[] = {NF, Nv, endBandSlice};

    auto CGajSliced(NEW(Tensor<complex>, sCGai->slice(CGajStart, CGajEnd)));
    auto conjCGajSliced(
        NEW(Tensor<complex>, conjCGai->slice(CGajStart, CGajEnd)));

    // Oldschool
    //    (*FGij)["GFij"] = (*conjCGai)["Gai"] * (*CGajSliced)["Faj"];
    //    (*FGji)["GFij"] = (*conjCGai)["Fbi"] * (*CGajSliced)["Gbj"];

    //    (*QGGs)["GF"]  += (0.5)* (*FGij)["GFij"] * (*FGji)["GFij"];
    //    (*QGGt)["GF"]  += (-0.75)* (*FGij)["GFij"] * (*FGji)["GFij"];
    // newschool
    (*FGij)["GFij"] = (*CGajSliced)["Gai"] * (*CGajSliced)["Faj"];
    (*FGji)["GFij"] = (*conjCGajSliced)["Gbj"] * (*conjCGajSliced)["Fbi"];
    (*QGGs)["GF"] += (0.5) * (*FGij)["GFij"] * (*FGji)["GFij"];
    (*QGGt)["GF"] += (-0.75) * (*FGij)["GFij"] * (*FGji)["GFij"];
    if (out.present("QGGx")) {
      (*cQGGx)["GF"] = (0.5) * (*FGij)["GFij"] * (*FGji)["GFij"];
    }
  }

  auto realQGGs(new Tensor<double>(2, NGG));
  auto realQGGt(new Tensor<double>(2, NGG));

  fromComplexTensor(*QGGs, *realQGGs);
  fromComplexTensor(*QGGt, *realQGGt);
  out.set<Tensor<double> *>("QGGs", realQGGs);
  out.set<Tensor<double> *>("QGGt", realQGGt);

  if (out.present("QGGx")) {
    auto QGGx(new Tensor<double>(2, NGG));
    fromComplexTensor(*cQGGx, *QGGx);
    out.set<Tensor<double> *>("QGGx", QGGx);
  }
}

void BasisSetExtrapolation::calculateNewSF(int type,
                                           real gamma,
                                           Tensor<double> *coulombKernel,
                                           Tensor<double> *newSF,
                                           Tensor<double> *resNewSF) {

  Tensor<double> *QGG(in.get<Tensor<double> *>("QGG"));
  int NG(QGG->lens[0]);
  int NFF[] = {NG};

  auto cK(new Tensor<double>(1, NFF));
  (*cK) = (*coulombKernel);

  real volume(in.get<double>("volume"));

  auto reciprocalYC(new Tensor<double>(1, NFF));

  CTF::Transform<real, real>(std::function<void(real, real &)>(
      [&type, &volume, &gamma](real cK, real &YC) {
        if (cK == 0) {
          YC = 0.;
        } else {
          if (type == 1) {
            YC = 1. / cK
               / cK; //*volume/4.5835494674469; //Ha(eV)*a0(A)*4*pi/(2*pi)**2
                     //        rYC = rYC/150.4121/(gamma*gamma+150.4121/rYC);
            //        hbar*hbar/2*me*10**20meter*(2pi)**2 in eV rYC =
            //        rYC*4.5835494674469/volume/(gamma*gamma+2.*150.4121/rYC);
            //        // GAMMA[Energy]
            YC *=
                (-0.015237) / volume / (gamma * gamma + 1. / YC); // GAMMA[1/A]
          } else if (type == 2) {
            YC = (-0.015237) / volume / (gamma * gamma + cK * cK)
               / (gamma * gamma + cK * cK);
          }
        }
      }))((*cK)["G"], (*reciprocalYC)["G"]);

  auto dReciprocalYC(new Tensor<double>(1, NFF));

  CTF::Transform<real, real, real>(std::function<void(real, real, real &)>(
      [&type, &volume, &gamma](real cK, real YC, real &dYC) {
        if (cK == 0) {
          dYC = 0.;
        } else {
          if (type == 1) {
            dYC = (-2.) * YC * gamma / (gamma * gamma + cK * cK);
          } else if (type == 2) {
            dYC = (-4.) * YC * gamma / (gamma * gamma + cK * cK);
          }
        }
      }))((*cK)["G"], (*reciprocalYC)["G"], (*dReciprocalYC)["G"]);

  (*newSF)["G"] = (*QGG)["GF"] * (*reciprocalYC)["F"];
  (*resNewSF)["G"] = (*QGG)["GF"] * (*dReciprocalYC)["F"];
}

void BasisSetExtrapolation::fitF12(int type, real minG, real maxG) {

  real volume(in.get<double>("volume"));
  if (volume < 0.) throw new EXCEPTION("Set volume");

  Tensor<double> *structureFactor(in.get<Tensor<double> *>("StructureFactor"));
  Tensor<double> *coulombKernel(in.get<Tensor<double> *>("CoulombKernel"));

  int NG(coulombKernel->lens[0]);
  int NFF[] = {NG};

  auto residuumFittedSF(new Tensor<double>(1, NFF));
  auto fittedSF(new Tensor<double>(1, NFF));
  auto absoluteG(new Tensor<double>(1, NFF));

  // Take out infinity
  CTF::Transform<real>(std::function<void(real &)>([](real &cK) {
    if (std::isinf(cK)) { cK = 0.; }
  }))((*coulombKernel)["G"]);
  // Construct array: aboluteG
  CTF::Transform<real, real>(
      std::function<void(real, real &)>([&volume](real cK, real &absG) {
        if (cK == 0.) {
          absG = 0.;
        } else {
          absG = cK * volume / 4.5835494674469;
          absG = 1. / std::sqrt(absG);
        }
      }))((*coulombKernel)["G"], (*absoluteG)["G"]);
  real gamma(in.get<double>("gamma"));
  int iterations(in.get<int64_t>("iterations"));
  for (int i(0); i <= iterations; ++i) {

    calculateNewSF(type, gamma, absoluteG, fittedSF, residuumFittedSF);

    auto dummy(new Tensor<double>(1, NFF));
    CTF::Transform<real, real, real>(std::function<void(real, real, real &)>(
        [&maxG, &minG](real absG, real res, real &dummy) {
          if (absG > maxG || absG < minG) {
            dummy = 0.;
          } else {
            dummy = res * res;
          }
        }))((*absoluteG)["G"], (*residuumFittedSF)["G"], (*dummy)["G"]);

    CTF::Scalar<double> denom;
    denom[""] = (*dummy)["G"];

    CTF::Transform<real, real, real>(std::function<void(real, real, real &)>(
        [](real fitsf, real sf, real &dummy) {
          dummy = fitsf - sf;
        }))((*fittedSF)["G"], (*structureFactor)["G"], (*dummy)["G"]);

    CTF::Transform<real, real, real>(std::function<void(real, real, real &)>(
        [&maxG, &minG](real absG, real resfitsf, real &dummy) {
          if (absG > maxG || absG < minG) {
            dummy = 0.;
          } else {
            dummy *= resfitsf;
          }
        }))((*absoluteG)["G"], (*residuumFittedSF)["G"], (*dummy)["G"]);

    CTF::Scalar<double> numer;
    numer[""] = (*dummy)["G"];
    gamma -= numer.get_val() / denom.get_val();

    (*residuumFittedSF) = (*structureFactor);

    CTF::Transform<real, real, real>(std::function<void(real, real, real &)>(
        [&maxG, &minG](real absG, real fitsf, real &resfitsf) {
          if (absG > maxG || absG < minG) {
            resfitsf = 0.;
          } else {
            resfitsf -= fitsf;
          }
        }))((*absoluteG)["G"], (*fittedSF)["G"], (*residuumFittedSF)["G"]);
    double residuum(residuumFittedSF->norm2());
    LOG(0, "gamma") << gamma << " norm: " << residuum << std::endl;
  }
  //  out.set<Tensor<double> *>("ResNewSF",resNewSF);
  out.set<Tensor<double> *>("FittedSF", fittedSF);
  // evaluate energy correction term E = v(G)*Q(G,G')*f12(G')
  CTF::Scalar<double> f12EnergyCorrection(*Sisi4s::world);
  f12EnergyCorrection[""] = (*coulombKernel)["G"] * (*fittedSF)["G"];
  double correction = f12EnergyCorrection.get_val();
  LOG(0, "gamma") << "f12 energy correction: " << correction << std::endl;
}

void BasisSetExtrapolation::invertQGG() {

  Tensor<double> *fullQGG(in.get<Tensor<double> *>("QGG"));

  int NG(fullQGG->lens[0]);
  int twodStart[] = {1, 1};
  int twodEnd[] = {NG, NG};
  Tensor<double> QGG(fullQGG->slice(twodStart, twodEnd));

  auto invQGG(new Tensor<double>(false, QGG));
  //  (*invQGG)["PQ"] = IterativePseudoInverse<complex>(QGG).get()["PQ"];
  (*invQGG)["PQ"] = PseudoInverseSvd<double>(QGG).get()["PQ"];
  int NF[] = {NG};
  Tensor<double> getStructureFactor(
      in.get<Tensor<double> *>("StructureFactor"));
  int NS(getStructureFactor.lens[0]);
  auto structureFactor(new Tensor<double>(1, NF));
  if (NS == NG) {
    LOG(0, "length") << NS << " " << NG << std::endl;
    (*structureFactor)["G"] = getStructureFactor["G"];
  } else if (NS == (NG + 1) / 2) {
    // structureFactor is halfmesh. however Q(G,G') is full mesh
    LOG(0, "length") << NS << " " << NG << std::endl;
    int dstStart[] = {0};
    int dstEnd[] = {NS};
    int srcStart[] = {0};
    int srcEnd[] = {NS};
    structureFactor->slice(dstStart,
                           dstEnd,
                           1.0,
                           getStructureFactor,
                           srcStart,
                           srcEnd,
                           0.5);
    dstStart[0] = NS;
    dstEnd[0] = NG;
    srcStart[0] = 1;
    srcEnd[0] = NS;
    structureFactor->slice(dstStart,
                           dstEnd,
                           1.0,
                           getStructureFactor,
                           srcStart,
                           srcEnd,
                           0.5);
  } else {
    LOG(0, "length") << NS << " " << NG << std::endl;
    throw new EXCEPTION("dimension problems of Q(G,G') and S(G)");
  }

  int onedStart[] = {1};
  int onedEnd[] = {NG};
  int NGG(NG - 1);
  int NFF[] = {NGG};

  auto slicedStructureFactor(
      new Tensor<double>(structureFactor->slice(onedStart, onedEnd)));

  auto slicedF12(new Tensor<double>(false, *slicedStructureFactor));

  (*slicedF12)["P"] = (*invQGG)["PQ"] * (*slicedStructureFactor)["Q"];

  auto f12(new Tensor<double>(1, NF));
  int dstStart[] = {0};
  f12->slice(onedStart, onedEnd, 1.0, *slicedF12, dstStart, NFF, 1.0);

  out.set<Tensor<double> *>("f12", f12);
}
