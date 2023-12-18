#include <array>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <algorithms/Mp2NaturalOrbitals.hpp>
// #include <math/MathFunctions.hpp>
// #include <math/ComplexTensor.hpp>
// #include <DryTensor.hpp>
#include <extern/Lapack.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
// #include <Options.hpp>

using namespace sisi4s;

IMPLEMENT_EMPTY_DRYRUN(Mp2NaturalOrbitals) {}

// We follow Taube & Bartlett:
// doi:10.1135/cccc20050837
// Diagonalize single particle density matrix (1)
// truncate using threshold (2)
// Rotate Fock matrix (4)
// diagonalze truncated Fock matrix (5)
// obtain new coefficients

DEFSPEC(
    Mp2NaturalOrbitals,
    SPEC_IN({"occupationThreshold", SPEC_VALUE_DEF("TODO: DOC", double, 1e-16)},
            {"FnoAlpha", SPEC_VALUE("TODO: DOC", int64_t)},
            {"FnoBeta", SPEC_VALUE("TODO: DOC", int64_t)},
            {"FnoNumber", SPEC_VALUE("TODO: DOC", int64_t)},
            {"unrestricted", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
            {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"OrbitalCoefficients", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"ParticleEigenEnergies",
             SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"PPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<double> *)},
            {"Spins", SPEC_VARIN("TODO: DOC", Tensor<double> *)}),
    SPEC_OUT({"occupationNumber", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"ParticleEigenEnergiesRediag",
              SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"ParticleEigenEnergiesRediag",
              SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"RotatedOrbitals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)},
             {"RotatedOrbitals", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(Mp2NaturalOrbitals) {
  Tensor<double> *orbs(in.get<Tensor<double> *>("OrbitalCoefficients"));
  Tensor<double> *Vabij(in.get<Tensor<double> *>("PPHHCoulombIntegrals"));
  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));
  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));
  const bool unrestricted(in.get<bool>("unrestricted"));
  auto rotatedOrbitals(new Tensor<double>(false, *orbs));
  auto epsaRediag(new Tensor<double>(false, *epsa));
  auto occNumber(new Tensor<double>(false, *epsa));
  int64_t Nv(epsa->lens[0]);
  int64_t No(epsi->lens[0]);
  int64_t Np(orbs->lens[1]);
  int64_t Nx(orbs->lens[0]);
  int64_t nFno;
  std::array<int, 2> vv({{(int)Nv, (int)Nv}});
  std::array<int, 2> syms({{NS, NS}});
  std::array<int, 2> pp({{(int)Np, (int)Np}});
  auto rotationMatrix(new Tensor<double>(2,
                                         pp.data(),
                                         syms.data(),
                                         *Vabij->wrld,
                                         "rotationMatrix"));
  auto Tabij(new Tensor<double>(Vabij));

  if (!unrestricted) {

    auto Tcbij(new Tensor<double>(Vabij));
    auto Dab(
        new Tensor<double>(2, vv.data(), syms.data(), *Vabij->wrld, "Dab"));

    Tabij->set_name("Tabij");
    Tcbij->set_name("Tcbij");
    (*Tabij)["abij"] = (*epsi)["i"];
    (*Tabij)["abij"] += (*epsi)["j"];
    (*Tabij)["abij"] -= (*epsa)["a"];
    (*Tabij)["abij"] -= (*epsa)["b"];
    CTF::Transform<double, double>(
        std::function<void(double, double &)>([](double vabij, double &tabij) {
          tabij = vabij / tabij;
        }))((*Vabij)["abij"], (*Tabij)["abij"]);
    (*Tcbij)["cbij"] = (2.0) * (*Tabij)["cbij"];
    (*Tcbij)["cbij"] += (-1.0) * (*Tabij)["cbji"];
    (*Dab)["ab"] = (*Tcbij)["caij"] * (*Tabij)["cbij"];

    std::vector<double> DabMatrix(Nv * Nv);
    std::vector<double> FabMatrix(Nv * Nv);
    std::vector<double> RabMatrix(Nv * Nv);
    std::vector<double> ParticleEigenEnergies(Nv);
    epsa->read_all(ParticleEigenEnergies.data());
    Dab->read_all(DabMatrix.data());

    int NvInt(Nv);
    int info;
    int lwork(-1);
    std::vector<double> w(Nv);
    double wlength;

    // DIAGONALIZE (1)
    dsyev_("V",
           "L",
           &NvInt,
           DabMatrix.data(),
           &NvInt,
           w.data(),
           &wlength,
           &lwork,
           &info);
    lwork = wlength;
    std::vector<double> work(lwork);
    dsyev_("V",
           "L",
           &NvInt,
           DabMatrix.data(),
           &NvInt,
           w.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization (1), naturalOrbitals\n";

    // TRUNCATE (2)
    double occupationThreshold(in.get<double>("occupationThreshold"));
    if (in.present("FnoNumber")) {
      int64_t fnoNumber(in.get<int64_t>("FnoNumber"));
      nFno = std::min(fnoNumber, Nv);
      // because of the 'wrong' ordering of the eigenvalues
      // we have to zero the first Nv-NFno columns
      for (int64_t a(0); a < Nv - fnoNumber; a++)
        for (int64_t b(0); b < Nv; b++) DabMatrix[b + a * Nv] = 0.0;
      LOG(0, "Nocc") << nFno << std::endl;
    } else {
      int counter(0);
      for (int64_t a(0); a < Nv; a++) {
        if (w[a] < occupationThreshold) {
          counter++;
          for (int64_t b(0); b < Nv; b++) DabMatrix[b + a * Nv] = 0.0;
        }
      }
      nFno = Nv - counter;
      LOG(0, "occupationThreshold") << occupationThreshold << std::endl;
      LOG(0, "Nfno") << Nv - counter << std::endl;
    }

    for (int64_t a(0); a < Nv; a++) LOG(2, "occ") << w[a] << std::endl;

    const int rank_m = int(Sisi4s::world->rank == 0); // rank mask
    std::vector<int64_t> index;
    index.resize(rank_m * Nv);
    std::iota(index.begin(), index.end(), 0);
    if (in.present("occupationNumber")) {
      occNumber->write(index.size(), index.data(), w.data());
      out.set<Tensor<double> *>("occupationNumber", occNumber);
      LOG(0, "writing:") << "occupationNumber\n";
    }

    // ROTATION (4)
    for (int64_t b(0); b < Nv; b++)
      for (int64_t a(0); a < Nv; a++)
        for (int64_t c(0); c < Nv; c++) {
          FabMatrix[a + b * Nv] += DabMatrix[c + a * Nv] * DabMatrix[c + b * Nv]
                                 * ParticleEigenEnergies[c];
        }

    dsyev_("V",
           "L",
           &NvInt,
           FabMatrix.data(),
           &NvInt,
           w.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization (4), naturalOrbitals\n";

    for (int64_t a(0); a < Nv; a++) LOG(2, "eVal") << w[a] << std::endl;

    // rotation matrix for orbtial (virtual) coefficients Dab*Fab
    for (int64_t b(0); b < Nv; b++)
      for (int64_t a(0); a < Nv; a++)
        for (int64_t c(0); c < Nv; c++) {
          RabMatrix[a + b * Nv] +=
              DabMatrix[a + c * Nv] * FabMatrix[c + b * Nv];
        }

    epsaRediag->write(index.size(), index.data(), w.data());

    std::vector<double> unity(No * No);
    for (int64_t i(0); i < No; i++) unity[i + i * No] = 1.0;

    std::array<int, 2> oo({{(int)No, (int)No}});
    auto newunity(
        new Tensor<double>(2, oo.data(), syms.data(), *Vabij->wrld, "new"));
    index.resize(rank_m * No * No);
    std::iota(index.begin(), index.end(), 0);
    newunity->write(index.size(), index.data(), unity.data());

    Tensor<double> virtualRotor(false, *Dab);

    index.resize(rank_m * Nv * Nv);
    std::iota(index.begin(), index.end(), 0);
    virtualRotor.write(index.size(), index.data(), RabMatrix.data());

    int dstStart[] = {0, 0};
    int dstEnd[] = {(int)No, (int)No};
    int srcStart[] = {0, 0};
    int srcEnd[] = {(int)No, (int)No};
    rotationMatrix
        ->slice(dstStart, dstEnd, 1.0, newunity, srcStart, srcEnd, 1.0);
    srcEnd[0] = Nv;
    srcEnd[1] = Nv;
    dstStart[0] = No;
    dstStart[1] = No;
    dstEnd[0] = No + Nv;
    dstEnd[1] = No + Nv;
    rotationMatrix
        ->slice(dstStart, dstEnd, 1.0, virtualRotor, srcStart, srcEnd, 1.0);
    LOG(0, "dims") << rotatedOrbitals->lens[0] << "x"
                   << rotatedOrbitals->lens[1] << "  =  " << orbs->lens[0]
                   << "x" << orbs->lens[1] << "  *  " << rotationMatrix->lens[0]
                   << "x" << rotationMatrix->lens[1] << std::endl;
    (*rotatedOrbitals)["mi"] = (*orbs)["mj"] * (*rotationMatrix)["ji"];

    std::array<int, 2> ff(
        {{static_cast<int>(Nx), static_cast<int>(No + nFno)}});
    std::array<int, 2> f({{static_cast<int>(nFno)}});

    auto newEpsa(
        new Tensor<double>(1, f.data(), syms.data(), *Vabij->wrld, "newEpsa"));
    dstStart[0] = 0;
    dstEnd[0] = nFno;
    srcStart[0] = Nv - nFno;
    srcEnd[1] = Nv;
    newEpsa->slice(dstStart, dstEnd, 1.0, epsaRediag, srcStart, srcEnd, 1.0);

    auto newRotatedOrbitals(
        new Tensor<double>(2, ff.data(), syms.data(), *Vabij->wrld, "R"));
    srcStart[0] = 0;
    srcStart[1] = 0;
    srcEnd[0] = Nx;
    srcEnd[1] = No;
    dstStart[0] = 0;
    dstStart[1] = 0;
    dstEnd[0] = Nx;
    dstEnd[1] = No;
    newRotatedOrbitals
        ->slice(dstStart, dstEnd, 1.0, rotatedOrbitals, srcStart, srcEnd, 1.0);

    srcStart[1] = No + Nv - nFno;
    srcEnd[1] = No + Nv;
    dstStart[1] = No;
    dstEnd[1] = nFno + No;
    newRotatedOrbitals
        ->slice(dstStart, dstEnd, 1.0, rotatedOrbitals, srcStart, srcEnd, 1.0);

    out.set<Tensor<double> *>("RotatedOrbitals", newRotatedOrbitals);
    out.set<Tensor<double> *>("ParticleEigenEnergiesRediag", newEpsa);
  } else { // UNRESTRICTED

    auto Spins(in.get<Tensor<double> *>("Spins"));
    auto Dab(
        new Tensor<double>(2, vv.data(), syms.data(), *Vabij->wrld, "Dab"));

    Tabij->set_name("Tabij");
    (*Tabij)["abij"] = (*epsi)["i"];
    (*Tabij)["abij"] += (*epsi)["j"];
    (*Tabij)["abij"] -= (*epsa)["a"];
    (*Tabij)["abij"] -= (*epsa)["b"];
    CTF::Transform<double, double>(
        std::function<void(double, double &)>([](double vabij, double &tabij) {
          tabij = vabij / tabij;
        }))((*Vabij)["abij"], (*Tabij)["abij"]);
    (*Dab)["ab"] = (*Tabij)["caij"] * (*Tabij)["cbij"];

    std::vector<double> FabMatrix(Nv * Nv);
    std::vector<double> RabMatrix(Nv * Nv);
    std::vector<double> ParticleEigenEnergies(Nv);
    epsa->read_all(ParticleEigenEnergies.data());

    int Nalpha(0), Nbeta(0);
    std::vector<double> spins(Nv + No);
    Spins->read_all(spins.data());
    for (size_t ii(No); ii < spins.size(); ii++) {
      if (spins[ii] > 0.0) Nalpha++;
      else Nbeta++;
    }

    int aStart[] = {0, 0};
    int aEnd[] = {Nalpha, Nalpha};
    auto ctfDalpha(Dab->slice(aStart, aEnd));
    auto ctfDbeta(Dab->slice(aEnd, vv.data()));

    std::vector<double> Dalpha(Nalpha * Nalpha);
    std::vector<double> Dbeta(Nbeta * Nbeta);
    std::vector<double> Falpha(Nalpha * Nalpha);
    std::vector<double> Fbeta(Nbeta * Nbeta);
    std::vector<double> Ralpha(Nalpha * Nalpha);
    std::vector<double> Rbeta(Nbeta * Nbeta);

    std::vector<double> epsaAlpha(Nalpha);
    std::vector<double> epsaBeta(Nbeta);

    for (int ii(0); ii < Nalpha; ii++)
      epsaAlpha[ii] = ParticleEigenEnergies[ii];

    for (int ii(0); ii < Nbeta; ii++)
      epsaBeta[ii] = ParticleEigenEnergies[ii + Nalpha];

    ctfDalpha.read_all(Dalpha.data());
    ctfDbeta.read_all(Dbeta.data());

    int info;
    int lwork(-1);
    std::vector<double> walpha(Nalpha);
    std::vector<double> wbeta(Nbeta);
    double wlength;

    // DIAGONALIZE (1)
    dsyev_("V",
           "L",
           &Nalpha,
           Dalpha.data(),
           &Nalpha,
           walpha.data(),
           &wlength,
           &lwork,
           &info);
    lwork = wlength;
    std::vector<double> work(lwork);
    dsyev_("V",
           "L",
           &Nalpha,
           Dalpha.data(),
           &Nalpha,
           walpha.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization (1) unatrualOrbitals\n";

    // TRUNCATE (2)
    double occupationThreshold(in.get<double>("occupationThreshold"));
    if (in.present("FnoAlpha")) {
      int fnoNumber(in.get<int64_t>("FnoAlpha"));
      // because of the 'wrong' ordering of the eigenvalues
      // we have to zero the first Nv-NFno columns
      for (int64_t a(0); a < Nalpha - fnoNumber; a++)
        for (int64_t b(0); b < Nalpha; b++) Dalpha[b + a * Nalpha] = 0.0;

      LOG(0, "NFnoAlpha") << fnoNumber << std::endl;
    } else {
      int counter(0);
      for (int64_t a(0); a < Nalpha; a++) {
        if (walpha[a] < occupationThreshold) {
          counter++;
          for (int64_t b(0); b < Nalpha; b++) Dalpha[b + a * Nalpha] = 0.0;
        }
      }
      LOG(0, "occupationThreshold") << occupationThreshold << std::endl;
      LOG(0, "NFnoAlpha") << Nalpha - counter << std::endl;
    }

    for (int64_t a(0); a < Nalpha; a++)
      LOG(2, "Alpha occ") << walpha[a] << std::endl;

    lwork = -1;

    ////////////////////
    /// BETA CHANNEL
    ////////////////////

    // DIAGONALIZE (1)
    dsyev_("V",
           "L",
           &Nbeta,
           Dbeta.data(),
           &Nbeta,
           wbeta.data(),
           &wlength,
           &lwork,
           &info);
    lwork = wlength;
    work.resize(lwork);
    dsyev_("V",
           "L",
           &Nbeta,
           Dbeta.data(),
           &Nbeta,
           wbeta.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization (2), unatural orbitals\n";

    // TRUNCATE (2)
    if (in.present("FnoBeta")) {
      int fnoNumber(in.get<int64_t>("FnoBeta"));
      // because of the 'wrong' ordering of the eigenvalues
      // we have to zero the first Nv-NFno columns
      for (int64_t a(0); a < Nbeta - fnoNumber; a++)
        for (int64_t b(0); b < Nbeta; b++) Dbeta[b + a * Nbeta] = 0.0;

      LOG(0, "NFnoBeta") << fnoNumber << std::endl;
    } else {
      int counter(0);
      for (int64_t a(0); a < Nbeta; a++) {
        if (wbeta[a] < occupationThreshold) {
          counter++;
          for (int64_t b(0); b < Nv; b++) Dbeta[b + a * Nbeta] = 0.0;
        }
      }
      LOG(0, "occupationThreshold") << occupationThreshold << std::endl;
      LOG(0, "NFnoBeta") << Nbeta - counter << std::endl;
    }

    for (int64_t a(0); a < Nbeta; a++)
      LOG(2, "Beta occ") << wbeta[a] << std::endl;

    // ROTATION (4)
    for (int64_t b(0); b < Nalpha; b++)
      for (int64_t a(0); a < Nalpha; a++)
        for (int64_t c(0); c < Nalpha; c++) {
          Falpha[a + b * Nalpha] +=
              Dalpha[c + a * Nalpha] * Dalpha[c + b * Nalpha] * epsaAlpha[c];
        }

    lwork = -1;
    dsyev_("V",
           "L",
           &Nalpha,
           Falpha.data(),
           &Nalpha,
           walpha.data(),
           &wlength,
           &lwork,
           &info);
    lwork = wlength;
    work.resize(lwork);

    dsyev_("V",
           "L",
           &Nalpha,
           Falpha.data(),
           &Nalpha,
           walpha.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization, (4) unatural orbtials\n";

    for (int64_t a(0); a < Nalpha; a++)
      LOG(2, "Alpha eVal") << walpha[a] << std::endl;

    ////////////////
    // BETA CHANNEL
    ///////////////

    for (int64_t b(0); b < Nbeta; b++)
      for (int64_t a(0); a < Nbeta; a++)
        for (int64_t c(0); c < Nbeta; c++) {
          Fbeta[a + b * Nbeta] +=
              Dbeta[c + a * Nbeta] * Dbeta[c + b * Nbeta] * epsaBeta[c];
        }

    lwork = -1;
    dsyev_("V",
           "L",
           &Nbeta,
           Fbeta.data(),
           &Nbeta,
           wbeta.data(),
           &wlength,
           &lwork,
           &info);
    lwork = wlength;
    work.resize(lwork);

    dsyev_("V",
           "L",
           &Nbeta,
           Fbeta.data(),
           &Nbeta,
           wbeta.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0) throw "problem diagonalization, (4) unatural orbitals\n";

    for (int64_t a(0); a < Nbeta; a++)
      LOG(2, "Beta eVal") << wbeta[a] << std::endl;

    // rotation matrix for orbtial (virtual) coefficients Dab*Fab
    for (int64_t b(0); b < Nalpha; b++)
      for (int64_t a(0); a < Nalpha; a++)
        for (int64_t c(0); c < Nalpha; c++) {
          Ralpha[a + b * Nalpha] +=
              Dalpha[a + c * Nalpha] * Falpha[c + b * Nalpha];
        }

    // rotation matrix for orbtial (virtual) coefficients Dab*Fab
    for (int64_t b(0); b < Nbeta; b++)
      for (int64_t a(0); a < Nbeta; a++)
        for (int64_t c(0); c < Nbeta; c++) {
          Rbeta[a + b * Nbeta] += Dbeta[a + c * Nbeta] * Fbeta[c + b * Nbeta];
        }

    // construct combined ParticleEigenEnergies
    const int rank_m = int(Sisi4s::world->rank == 0); // rank mask
    std::vector<int64_t> index;
    index.resize(rank_m * Nv);
    std::iota(index.begin(), index.end(), 0);

    std::vector<double> w(Nv);
    for (int ii(0); ii < Nalpha; ii++) w[ii] = walpha[ii];
    for (int ii(0); ii < Nbeta; ii++) w[ii + Nalpha] = wbeta[ii];

    epsaRediag->write(index.size(), index.data(), w.data());

    // construct No-No block of transformation matrix
    std::vector<double> unity(No * No);
    for (int64_t i(0); i < No; i++) unity[i + i * No] = 1.0;

    std::array<int, 2> oo({{(int)No, (int)No}});
    auto newunity(
        new Tensor<double>(2, oo.data(), syms.data(), *Vabij->wrld, "new"));
    index.resize(rank_m * No * No);
    std::iota(index.begin(), index.end(), 0);
    newunity->write(index.size(), index.data(), unity.data());

    std::array<int, 2> aa({{Nalpha, Nalpha}});
    auto alphaRotor(
        new Tensor<double>(2, aa.data(), syms.data(), *Vabij->wrld, "aRotor"));
    index.resize(rank_m * Nalpha * Nalpha);
    std::iota(index.begin(), index.end(), 0);
    alphaRotor->write(index.size(), index.data(), Ralpha.data());

    std::array<int, 2> bb({{Nbeta, Nbeta}});
    auto betaRotor(
        new Tensor<double>(2, bb.data(), syms.data(), *Vabij->wrld, "aRotor"));
    index.resize(rank_m * Nbeta * Nbeta);
    std::iota(index.begin(), index.end(), 0);
    betaRotor->write(index.size(), index.data(), Rbeta.data());

    int dstStart[] = {0, 0};
    int dstEnd[] = {(int)No, (int)No};
    int srcStart[] = {0, 0};
    int srcEnd[] = {(int)No, (int)No};
    rotationMatrix
        ->slice(dstStart, dstEnd, 1.0, newunity, srcStart, srcEnd, 1.0);
    srcEnd[0] = Nalpha;
    srcEnd[1] = Nalpha;
    dstStart[0] = No;
    dstStart[1] = No;
    dstEnd[0] = No + Nalpha;
    dstEnd[1] = No + Nalpha;
    rotationMatrix
        ->slice(dstStart, dstEnd, 1.0, alphaRotor, srcStart, srcEnd, 1.0);
    srcEnd[0] = Nbeta;
    srcEnd[1] = Nbeta;
    dstStart[0] = No + Nalpha;
    dstStart[1] = No + Nalpha;
    dstEnd[0] = Np;
    dstEnd[1] = Np;
    rotationMatrix
        ->slice(dstStart, dstEnd, 1.0, betaRotor, srcStart, srcEnd, 1.0);

    LOG(0, "dims") << rotatedOrbitals->lens[0] << "x"
                   << rotatedOrbitals->lens[1] << "  =  " << orbs->lens[0]
                   << "x" << orbs->lens[1] << "  *  " << rotationMatrix->lens[0]
                   << "x" << rotationMatrix->lens[1] << std::endl;
    (*rotatedOrbitals)["mi"] = (*orbs)["mj"] * (*rotationMatrix)["ji"];

    out.set<Tensor<double> *>("RotatedOrbitals", rotatedOrbitals);
    out.set<Tensor<double> *>("ParticleEigenEnergiesRediag", epsaRediag);
  }
}
