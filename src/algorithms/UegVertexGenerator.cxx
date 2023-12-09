#include <algorithms/UegVertexGenerator.hpp>

#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <math/Complex.hpp>

using namespace sisi4s;

static double evalMadelung(const double v) {
  double kappa = pow(v, -1.0 / 3.0);
  double term2 = M_PI / (kappa * kappa * v);
  double term4 = 2 * kappa / sqrt(M_PI);
  double boxLength = 1.0 / kappa;
  double recipsum = 0.0;
  double realsum = 0.0;
  for (int l1 = -6; l1 <= 6; ++l1)
    for (int l2 = -6; l2 <= 6; ++l2)
      for (int l3 = -6; l3 <= 6; ++l3) {
        int n2 = l1 * l1 + l2 * l2 + l3 * l3;
        double modr = boxLength * sqrt((double)n2);
        double k2 = kappa * kappa * n2;
        if (n2 > 0) {
          recipsum -=
              1.0 / (M_PI * k2) * exp(-M_PI * M_PI * k2 / kappa / kappa) / v;
          realsum -= erfc(kappa * modr) / modr;
        }
      }
  return realsum + term2 + term4 + recipsum;
}

// define two functions which give the squared length of the grid-points
static size_t sL(const ivec a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}
static double sL(const dvec a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}
static double
Vijji(const double madelung, const dvec a, const dvec b, const double v) {
  dvec q({a[0] - b[0], a[1] - b[1], a[2] - b[2]});
  if (sL(q) < 1e-8) return madelung;
  return 4.0 * M_PI / v / sL(q);
}

IMPLEMENT_EMPTY_DRYRUN(UegVertexGenerator) {}

DEFSPEC(
    UegVertexGenerator,
    SPEC_IN(
        {"madelung",
         SPEC_VALUE_DEF("The approximation for V(G=0)", double, -1.0)},
        {"rs", SPEC_POSITIVE("Electron density to be used", double)},
        {"halfGrid", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
        {"hartreeFock",
         SPEC_VALUE_DEF("Wether or not to calculate Hartree-Fock approximation",
                        bool,
                        true)},
        {"NF", SPEC_VALUE_DEF("TODO: DOC", size_t, 0)},
        {"No", SPEC_POSITIVE("The number of occupied orbitals", size_t)},
        {"Nv", SPEC_VALUE("The number of unoccupied orbitals", size_t)}),
    SPEC_OUT({"HoleEigenEnergies",
              SPEC_VAROUT("The obtained one-body hole energies",
                          Tensor<double> *)},
             {"ParticleEigenEnergies",
              SPEC_VAROUT("The obtained one-body particle energies",
                          Tensor<double> *)},
             {"CoulombVertex",
              SPEC_VAROUT("The Coulomb vergex Γ[Gpq] calculated",
                          Tensor<sisi4s::complex> *)}));

IMPLEMENT_ALGORITHM(UegVertexGenerator) {
  // We use the HF reference by default.
  bool lhfref(in.get<bool>("hartreeFock"));
  bool lclosed(true);
  const size_t No = in.get<size_t>("No"), Nv = in.get<size_t>("Nv");
  size_t NF = in.get<size_t>("NF");
  double rs = in.get<double>("rs");
  halfGrid = in.get<bool>("halfGrid");
  double madelung = in.get<double>("madelung");
  size_t Np(No + Nv);
  if (rs <= 0.0) throw("Invalid rs");

  // setup the integer Grid.
  //  1) gather more than enough candidates
  //  2.) sort by length
  //  3.) split and cut
  int maxG = pow(5.0 * Np, 1.0 / 3.0);
  std::vector<ivec> iGrid;
  for (int g1(-maxG); g1 <= maxG; g1++)
    for (int g2(-maxG); g2 <= maxG; g2++)
      for (int g3(-maxG); g3 <= maxG; g3++) iGrid.push_back({g1, g2, g3});

  sort(iGrid.begin(), iGrid.end(), [](ivec a, ivec b) {
    return sL(a) < sL(b);
  });
  if (iGrid.size() < Np) { throw("BUG related to Np & maxG\n"); }
  if (sL(iGrid[No]) == sL(iGrid[No - 1])) {
    OUT() << "WARNING: occupied orbitals form not a closed shell\n";
    if (!lhfref)
      throw("Zero gap system! Either change No or use Hartree-Fock!");
    lclosed = false;
  }

  if (sL(iGrid[Np]) == sL(iGrid[Np - 1])) {
    OUT() << "WARNING: virtual orbitals form not a closed shell\n";
  }
  iGrid.resize(Np);

  // define volume, lattice Constant, and reciprocal lattice constant
  double v(rs * rs * rs / 3.0 * 4.0 * M_PI * No * 2);
  double a(pow(v, 1. / 3.));
  double b(2.0 * M_PI / a);

  if (madelung < 0.0) madelung = evalMadelung(v);

  std::vector<dvec> dGrid;
  // here we can introduce a possible shift of the mesh
  for (auto i : iGrid) dGrid.push_back({b * i[0], b * i[1], b * i[2], 0.0});

  // now we can write the hartree fock energy in the 4th entry
  for (auto &d : dGrid) {
    d[3] = 0.5 * sL(d); // add the kinetic energy
    double exchE(0.0);
    for (size_t o(0); o < No; o++) exchE += Vijji(madelung, d, dGrid[o], v);
    if (lhfref) d[3] -= exchE;
  }
  double refE(0.0);
  for (size_t o(0); o < No; o++) {
    refE += dGrid[o][3];
    if (lhfref) refE += 0.5 * sL(dGrid[o]);
  }

  std::vector<double> energies(dGrid.size());

  for (size_t d(0); d < dGrid.size(); d++) energies[d] = dGrid[d][3];

  // double fermiEnergy((energies[No] + energies[No - 1]) / 2.0);

  // construct the momentum transition grid
  // 1.) get the largest momentum vector between two states p - q
  // 2.) construct a full grid with a largest grid vec. of this size

  ivec maxMom({0, 0, 0});
  for (size_t p(0); p < Np; p++)
    for (size_t q(0); q < Np; q++) {
      ivec d = {iGrid[p][0] - iGrid[q][0],
                iGrid[p][1] - iGrid[q][1],
                iGrid[p][2] - iGrid[q][2]};
      maxMom = max(maxMom, d, [](ivec a, ivec b) { return sL(a) < sL(b); });
    }

  size_t maxR = sL(maxMom);
  maxG = max({maxMom[0], maxMom[1], maxMom[2]},
             [](int a, int b) { return std::abs(a) < std::abs(b); });
  maxG = std::abs(maxG);
  std::map<ivec, size_t> momMap;
  size_t index(0);
  for (int g1(-maxG); g1 <= maxG; g1++)
    for (int g2(-maxG); g2 <= maxG; g2++)
      for (int g3(-maxG); g3 <= maxG; g3++) {
        ivec t({g1, g2, g3});
        if (sL(t) > maxR) continue;
        momMap[t] = index++;
      }
  if (NF == 0) NF = momMap.size();

  if (NF != momMap.size() || halfGrid || !lclosed) {
    OUT() << "WARNING: the Vertex will not be correct! Just for profiling!\n";
  }

  double fac(4.0 * M_PI / v);

  int syms[] = {NS, NS, NS, NS};
  const int coulombVertexLens[] = {(int)NF, (int)Np, (int)Np};
  auto coulombVertex = new Tensor<complex>(3,
                                           coulombVertexLens,
                                           syms,
                                           *Sisi4s::world,
                                           "CoulombVertex");

  OUT() << "System Information:\n";
  OUT() << std::setprecision(3) << "  rs " << rs << ", No " << No << ", Nv "
        << Nv << "\n";
  OUT() << std::setprecision(10) << "  Volume " << v << ", madelung "
        << madelung << "\n";
  OUT() << "  HOMO " << energies[No - 1] << ", LUMO " << energies[No] << "\n";
  OUT() << "  Reference Energy per Electron/total " << refE / No / 2 << "/"
        << refE << std::endl;

  // Prepare eigenEnergies
  const int o[] = {(int)No}, _v[] = {(int)Nv};
  auto epsi = new Tensor<double>(1, o, syms, *Sisi4s::world, "epsi"),
       epsa = new Tensor<double>(1, _v, syms, *Sisi4s::world, "epsa");

  out.set<Tensor<sisi4s::complex> *>("CoulombVertex", coulombVertex);
  out.set<Tensor<double> *>("HoleEigenEnergies", epsi);
  out.set<Tensor<double> *>("ParticleEigenEnergies", epsa);

  // if we are in a dryRun there is nothing more to do
  // if (Sisi4s::dryRun) return result;

  size_t np = Sisi4s::world->np;
  size_t rank = Sisi4s::world->rank;
  // only rank 0 writes the data to the tensor
  std::vector<size_t> idx;
  if (!Sisi4s::world->rank) {
    idx.resize(energies.size());
    std::iota(idx.begin(), idx.end(), 0);
  }

  epsi->write(No, (int64_t *)idx.data(), energies.data());
  epsa->write(Nv, (int64_t *)idx.data() + No, energies.data());

  // Writing CoulombVertex to buffer
  // We have to do it mpi-able...otherwise we will
  // not be able to write it to a ctf tensor
  // We slice the number of states for all the mpi processes
  size_t slices(Np / np);
  std::vector<size_t> slicePerRank(np);
  for (size_t r(0); r < np; r++) {
    size_t lslice(slices);
    for (size_t i(0); i < Np - slices * np; i++)
      if (r == i) { lslice++; }
    slicePerRank[r] = lslice;
  }
  slices = slicePerRank[rank];
  // allocate only a buffer of needed size
  std::vector<complex> out(NF * Np * slices, {0, 0});
  // determine begin and end of the rank's slices
  auto sbegin(std::accumulate(slicePerRank.begin(),
                              slicePerRank.begin() + rank,
                              0UL,
                              std::plus<size_t>()));

  for (size_t s(0); s < slices; s++)
    for (size_t q(0); q < Np; q++) {
      auto p(s + sbegin);
      ivec d = {iGrid[q][0] - iGrid[p][0],
                iGrid[q][1] - iGrid[p][1],
                iGrid[q][2] - iGrid[p][2]};
      // This is a hack!
      // If NF is chosen by the user we will not have an overflow
      size_t ii = momMap[d] % NF;
      double res;
      (sL(d)) ? res = fac / (b * b * sL(d)) : res = evalMadelung(v);
      out[ii + q * NF + s * NF * Np] = {sqrt(res), 0.0};
    }

  idx.resize(out.size());
  std::iota(idx.begin(), idx.end(), sbegin * Np * NF);
  coulombVertex->write(idx.size(), (int64_t *)idx.data(), out.data());
}
