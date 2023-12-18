#include <algorithms/CcsdtEnergyFromCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>
#include <array>

using namespace sisi4s;

using F = double;
DEFSPEC(CcsdtEnergyFromCoulombIntegrals,
        SPEC_IN(CLUSTER_SINGLES_DOUBLES_TRIPLES_INSPEC,
                {"factorsSliceSize", SPEC_VALUE_DEF("TODO: DOC", int64_t, -1)},
                {"factorsSliceFactor", SPEC_VALUE("TODO: DOC", double)},
                {"integralsSliceFactor", SPEC_VALUE("TODO: DOC", double)},
                {"HHHHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HHHPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HHPHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HHPPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HPHHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HPHPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HPPHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"HPPPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PHHHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PHHPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PHPPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PPHPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PPPHCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PPPPCoulombIntegrals",
                 SPEC_VARIN("DOC: TODO", Tensor<F> *)->require()},
                {"PHPHCoulombIntegrals",
                 SPEC_VARIN("TODO: DOC", Tensor<F> *)->require()}),
        SPEC_OUT(CLUSTER_SINGLES_DOUBLES_TRIPLES_OUTSPEC));

template <typename F>
PTR(FockVector<F>) CcsdtEnergyFromCoulombIntegrals::getResiduumTemplate(
    const int iterationStep,
    const PTR(const FockVector<F>) &amplitudes) {

  auto *epsi = in.get<Tensor<double> *>("HoleEigenEnergies"),
       *epsa = in.get<Tensor<double> *>("ParticleEigenEnergies");

  auto Vhhhh(in.get<Tensor<F> *>("HHHHCoulombIntegrals")),
      Vpppp(in.get<Tensor<F> *>("PPPPCoulombIntegrals")),
      Vhhhp(in.get<Tensor<F> *>("HHHPCoulombIntegrals")),
      Vhhpp(in.get<Tensor<F> *>("HHPPCoulombIntegrals")),
      Vhphh(in.get<Tensor<F> *>("HPHHCoulombIntegrals")),
      Vhphp(in.get<Tensor<F> *>("HPHPCoulombIntegrals")),
      Vhppp(in.get<Tensor<F> *>("HPPPCoulombIntegrals")),
      Vpphh(in.get<Tensor<F> *>("PPHHCoulombIntegrals")),
      Vpphp(in.get<Tensor<F> *>("PPHPCoulombIntegrals")),
      Vhpph(in.get<Tensor<F> *>("HPPHCoulombIntegrals")),
      Vphpp(in.get<Tensor<F> *>("PHPPCoulombIntegrals")),
      Vhhph(in.get<Tensor<F> *>("HHPHCoulombIntegrals")),
      Vppph(in.get<Tensor<F> *>("PPPHCoulombIntegrals")),
      Vphph(in.get<Tensor<F> *>("PHPHCoulombIntegrals")),
      Vphhp(in.get<Tensor<F> *>("PHHPCoulombIntegrals")),
      Vphhh(in.get<Tensor<F> *>("PHHHCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  /*
  Tensor<F>
    *Fab = new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab"),
    *Fij = new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij");


  { // build Fab and Fij
    const auto
      to_F = CTF::Transform<double, F>([](double eps, F &f) {f = eps;});
    to_F((*epsi)["i"],
         (*Fij)["ii"]);
    to_F((*epsa)["a"],
         (*Fab)["aa"]);
  }
*/

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  auto Tph(amplitudes->get(0));
  Tph->set_name("Tph");
  auto Tpphh(amplitudes->get(1));
  Tpphh->set_name("Tpphh");
  auto Tppphhh(amplitudes->get(2));
  Tppphhh->set_name("Tppphhh");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rph(residuum->get(0));
  Rph->set_name("Rph");
  auto Rpphh(residuum->get(1));
  Rpphh->set_name("Rpphh");
  auto Rppphhh(residuum->get(2));
  Rppphhh->set_name("Rppphhh");

  if ((iterationStep == 0) && !in.present("initialDoublesAmplitudes")) {
    LOG(1, getAbbreviation())
        << "Set initial Rpphh amplitudes to Vijab" << std::endl;
    (*Rpphh)["abij"] = (*Vhhpp)["ijab"];
    return residuum;
  }

  // CCSD equations

  std::array<int, 4> syms({{NS, NS, NS, NS}});
  std::array<int, 4> voov({{Nv, No, No, Nv}});
  std::array<int, 2> vo({{Nv, No}});
  std::array<int, 2> vv({{Nv, Nv}});
  std::array<int, 2> oo({{No, No}});
  std::array<int, 2> ov({{No, Nv}});

  //*************************************************************************
  //****************  T2 amplitude equations  *******************************
  //*************************************************************************

  LOG(1, getCapitalizedAbbreviation())
      << "Solving T2 Amplitude Equations" << std::endl;

  // Define intermediates
  Tensor<double> Lac(2, vv.data(), syms.data(), *Vpphh->wrld, "Lac");
  Tensor<double> Kac(2, vv.data(), syms.data(), *Vpphh->wrld, "Kac");
  Tensor<double> Lki(2, oo.data(), syms.data(), *Vpphh->wrld, "Lki");
  Tensor<double> Kki(2, oo.data(), syms.data(), *Vpphh->wrld, "Kki");
  Tensor<double> Kck(2, vo.data(), syms.data(), *Vpphh->wrld, "Kck");

  Tensor<double> Xklij(false, *Vhhhh);
  Tensor<double> Xabij(false, *Vhhhh);
  Tensor<double> Xakci(false, *Vphph);
  Tensor<double> Xabcd(false, *Vpppp);
  Tensor<double> Xakic(4, voov.data(), syms.data(), *Vpphh->wrld, "Xakic");
  Xklij.set_name("Xklij");
  Xabij.set_name("Xabij");
  Xakci.set_name("Xakci");
  Xabcd.set_name("Xabcd");
  Xakic.set_name("Xakic");
  // Build Kac
  Kac["ac"] = (-2.0) * (*Vpphh)["cdkl"] * (*Tpphh)["adkl"];
  Kac["ac"] += (1.0) * (*Vpphh)["dckl"] * (*Tpphh)["adkl"];
  Kac["ac"] += (-2.0) * (*Vpphh)["cdkl"] * (*Tph)["ak"] * (*Tph)["dl"];
  Kac["ac"] += (1.0) * (*Vpphh)["dckl"] * (*Tph)["ak"] * (*Tph)["dl"];

  // Build Lac
  Lac["ac"] = Kac["ac"];
  Lac["ac"] += (2.0) * (*Vppph)["cdak"] * (*Tph)["dk"];
  Lac["ac"] += (-1.0) * (*Vppph)["dcak"] * (*Tph)["dk"];

  // Build Kki
  Kki["ki"] = (2.0) * (*Vpphh)["cdkl"] * (*Tpphh)["cdil"];
  Kki["ki"] += (-1.0) * (*Vpphh)["dckl"] * (*Tpphh)["cdil"];
  Kki["ki"] += (2.0) * (*Vpphh)["cdkl"] * (*Tph)["ci"] * (*Tph)["dl"];
  Kki["ki"] += (-1.0) * (*Vpphh)["dckl"] * (*Tph)["ci"] * (*Tph)["dl"];

  // Build Lki
  Lki["ki"] = Kki["ki"];
  Lki["ki"] += (2.0) * (*Vhhhp)["klic"] * (*Tph)["cl"];
  Lki["ki"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["cl"];

  // Contract Lac with T2 Amplitudes
  (*Rpphh)["abij"] += (1.0) * Lac["ac"] * (*Tpphh)["cbij"];

  // Contract Lki with T2 Amplitudes
  (*Rpphh)["abij"] += (-1.0) * Lki["ki"] * (*Tpphh)["abkj"];

  // Contract Coulomb integrals with T2 amplitudes
  (*Rpphh)["abij"] += (1.0) * (*Vppph)["baci"] * (*Tph)["cj"];
  (*Rpphh)["abij"] += (-1.0) * (*Vphph)["bkci"] * (*Tph)["ak"] * (*Tph)["cj"];
  (*Rpphh)["abij"] += (-1.0) * (*Vhhhp)["jika"] * (*Tph)["bk"];
  (*Rpphh)["abij"] += (-1.0) * (*Vpphh)["acik"] * (*Tph)["cj"] * (*Tph)["bk"];

  // Build Xakic
  Xakic["akic"] = (*Vpphh)["acik"];
  Xakic["akic"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["al"];
  Xakic["akic"] += (1.0) * (*Vppph)["acdk"] * (*Tph)["di"];
  Xakic["akic"] += (-0.5) * (*Vpphh)["dclk"] * (*Tpphh)["dail"];
  Xakic["akic"] += (-1.0) * (*Vpphh)["dclk"] * (*Tph)["di"] * (*Tph)["al"];
  Xakic["akic"] += (1.0) * (*Vpphh)["dclk"] * (*Tpphh)["adil"];
  Xakic["akic"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["adil"];

  // Build Xakci
  Xakci["akci"] = (*Vphph)["akci"];
  Xakci["akci"] += (-1.0) * (*Vhhhp)["klic"] * (*Tph)["al"];
  Xakci["akci"] += (1.0) * (*Vppph)["adck"] * (*Tph)["di"];
  Xakci["akci"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["dail"];
  Xakci["akci"] += (-1.0) * (*Vpphh)["cdlk"] * (*Tph)["di"] * (*Tph)["al"];

  // Contract Xakic and Xakci intermediates with T2 amplitudes Tabij
  (*Rpphh)["abij"] += (2.0) * Xakic["akic"] * (*Tpphh)["cbkj"];
  (*Rpphh)["abij"] += (-1.0) * Xakic["akic"] * (*Tpphh)["bckj"];

  (*Rpphh)["abij"] += (-1.0) * Xakci["akci"] * (*Tpphh)["cbkj"];
  (*Rpphh)["abij"] += (-1.0) * Xakci["bkci"] * (*Tpphh)["ackj"];

  // Symmetrize Rabij by applying permutation operator
  // to save memory we use Xakci as intermediate for the permutation operator
  Xakci["aibj"] = (*Rpphh)["abij"];
  (*Rpphh)["abij"] += Xakci["bjai"];

  //////////////////////////////////////////////////////////////////////
  // Now add all terms to Rabij that do not need to be symmetrized with
  // the permutation operator
  //////////////////////////////////////////////////////////////////////

  // Rabij are the Tabij amplitudes for the next iteration and need to be build
  (*Rpphh)["abij"] += (*Vpphh)["abij"];

  // Build Xklij intermediate
  Xklij["klij"] = (*Vhhhh)["klij"];
  Xklij["klij"] += (*Vhhhp)["klic"] * (*Tph)["cj"];
  Xklij["klij"] += (*Vhhhp)["lkjc"] * (*Tph)["ci"];
  Xklij["klij"] += (*Vpphh)["cdkl"] * (*Tpphh)["cdij"];
  Xklij["klij"] += (*Vpphh)["cdkl"] * (*Tph)["ci"] * (*Tph)["dj"];

  // Contract Xklij with T2 Amplitudes
  (*Rpphh)["abij"] += Xklij["klij"] * (*Tpphh)["abkl"];

  // Contract Xklij with T1 Amplitudes
  (*Rpphh)["abij"] += Xklij["klij"] * (*Tph)["ak"] * (*Tph)["bl"];

  // Build Xabcd intermediate
  Xabcd["abcd"] = (1.0) * (*Vpppp)["abcd"];
  Xabcd["abcd"] += (-1.0) * (*Vppph)["cdak"] * (*Tph)["bk"];
  Xabcd["abcd"] += (-1.0) * (*Vppph)["dcbk"] * (*Tph)["ak"];

  // Contract Xabcd with T2 and T1 Amplitudes
  (*Rpphh)["abij"] += Xabcd["abcd"] * (*Tpphh)["cdij"];
  (*Rpphh)["abij"] += Xabcd["abcd"] * (*Tph)["ci"] * (*Tph)["dj"];

  //********************************************************************************
  //***********************  T1 amplitude equations
  //*******************************
  //********************************************************************************
  LOG(1, getCapitalizedAbbreviation())
      << "Solving T1 Amplitude Equations" << std::endl;

  // Contract Kac and Kki with T1 amplitudes
  (*Rph)["ai"] += (1.0) * Kac["ac"] * (*Tph)["ci"];
  (*Rph)["ai"] += (-1.0) * Kki["ki"] * (*Tph)["ak"];

  // Build Kck
  Kck["ck"] = (2.0) * (*Vpphh)["cdkl"] * (*Tph)["dl"];
  Kck["ck"] += (-1.0) * (*Vpphh)["cdlk"] * (*Tph)["dl"];

  // Contract all the rest terms with T1 and T2 amplitudes
  (*Rph)["ai"] += (2.0) * Kck["ck"] * (*Tpphh)["caki"];
  (*Rph)["ai"] += (-1.0) * Kck["ck"] * (*Tpphh)["caik"];
  // TODO: check whether it's +
  (*Rph)["ai"] += (1.0) * Kck["ck"] * (*Tph)["ci"] * (*Tph)["ak"];
  (*Rph)["ai"] += (2.0) * (*Vpphh)["acik"] * (*Tph)["ck"];
  (*Rph)["ai"] += (-1.0) * (*Vphph)["akci"] * (*Tph)["ck"];
  (*Rph)["ai"] += (2.0) * (*Vppph)["cdak"] * (*Tpphh)["cdik"];
  (*Rph)["ai"] += (-1.0) * (*Vppph)["dcak"] * (*Tpphh)["cdik"];
  (*Rph)["ai"] += (2.0) * (*Vppph)["cdak"] * (*Tph)["ci"] * (*Tph)["dk"];
  (*Rph)["ai"] += (-1.0) * (*Vppph)["dcak"] * (*Tph)["ci"] * (*Tph)["dk"];
  (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tpphh)["ackl"];
  (*Rph)["ai"] += (1.0) * (*Vhhhp)["lkic"] * (*Tpphh)["ackl"];
  (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tph)["ak"] * (*Tph)["cl"];
  (*Rph)["ai"] += (1.0) * (*Vhhhp)["lkic"] * (*Tph)["ak"] * (*Tph)["cl"];

  // T3 -> R1

  (*Rph)["ai"] += (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"];
  (*Rph)["ai"] += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"];
  (*Rph)["ai"] += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"];
  (*Rph)["ai"] += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"];

  // T3 -> R2
  Tensor<double> Wphpp(false, *Vphpp);
  Tensor<double> Whhhp(false, *Vhhhp);
  Tensor<double> Whhpp(false, *Vhhpp);
  Tensor<double> Xpphh(false, *Vpphh);
  Wphpp.set_name("Wphpp");
  Whhhp.set_name("Whhhp");
  Whhpp.set_name("Whhpp");
  Xpphh.set_name("Xpphh");

  // Pure T3->R2
  Wphpp["akcd"] = (2.0) * (*Vphpp)["akcd"];
  Wphpp["akcd"] += (-1.0) * (*Vphpp)["akdc"];
  Whhhp["klic"] = (2.0) * (*Vhhhp)["klic"];
  Whhhp["klic"] += (-1.0) * (*Vhhhp)["lkic"];
  Whhpp["klcd"] = (2.0) * (*Vhhpp)["klcd"];
  Whhpp["klcd"] += (-1.0) * (*Vhhpp)["lkcd"];

  Xpphh["abij"] = Wphpp["akcd"] * (*Tppphhh)["cbdijk"];
  Xpphh["abij"] += (-1.0) * (*Vphpp)["akcd"] * (*Tppphhh)["cdbijk"];
  Xpphh["abij"] += (-1.0) * Whhhp["klic"] * (*Tppphhh)["abckjl"];
  Xpphh["abij"] += (*Vhhhp)["klic"] * (*Tppphhh)["acbkjl"];
  Xpphh["abij"] += Xpphh["baji"];
  (*Rpphh)["abij"] += Xpphh["abij"];
  // T1+T3 -> R2
  Xpphh["abij"] = Whhpp["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"];
  Xpphh["abij"] += (-1.0) * Whhpp["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"];
  Xpphh["abij"] += (-1.0) * Whhpp["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"];
  Xpphh["abij"] += (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"];
  Xpphh["abij"] += (-1.0) * Whhpp["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"];
  Xpphh["abij"] += (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"];
  Xpphh["abij"] += Xpphh["baji"];
  (*Rpphh)["abij"] += Xpphh["abij"];

  // T1+T2+T3->R3

  Tensor<double> Xabie(false, *Vpphp);
  Tensor<double> Xamij(false, *Vphhh);
  Tensor<double> Xjkmn(false, *Vhhhh);
  Tensor<double> Xbcef(false, *Vpppp);
  Tensor<double> Xamie(false, *Vphhp);
  Tensor<double> Xamei(false, *Vphph);
  Tensor<double> Xp3h3(false, *Tppphhh);
  Tensor<double> Fphpp(false, *Vphpp);
  Tensor<double> Fphhh(false, *Vphhh);
  Tensor<double> Fphhp(false, *Vphhp);
  Tensor<double> Fphph(false, *Vphph);
  Tensor<double> Fhpph(false, *Vhpph);
  Tensor<double> Fhphp(false, *Vhphp);
  Tensor<double> Fhp(2, ov.data(), syms.data(), *Vpphh->wrld, "Fhp");
  Tensor<double> Xim(false, Lki);
  Tensor<double> Xae(false, Lac);

  Xabie.set_name("Xabie");
  Xamij.set_name("Xamij");
  Xjkmn.set_name("Xjkmn");
  Xbcef.set_name("Xbcef");
  Xamie.set_name("Xamie");
  Xamei.set_name("Xamei");
  Xp3h3.set_name("Xp3h3");
  Fphpp.set_name("Fphpp");
  Fphhh.set_name("Fphhh");
  Fphhp.set_name("Fphhp");
  Fphph.set_name("Fphph");
  Fhpph.set_name("Fhpph");
  Fhphp.set_name("Fhphp");
  Xim.set_name("Xim");
  Xae.set_name("Xae");

  // II.16
  Xpphh["abij"] = (*Tpphh)["abij"];
  Xpphh["abij"] += (*Tph)["ai"] * (*Tph)["bj"];

  // II.9
  Fphpp["amef"] = (*Vphpp)["amef"];
  Fphpp["amef"] += (-1.0) * (*Vhhpp)["nmef"] * (*Tph)["an"];

  // II.10
  Fphhh["eimn"] = (*Vphhh)["eimn"];
  Fphhh["eimn"] += (*Vhhpp)["mnef"] * (*Tph)["fi"];

  // II.11
  Fphhp["amie"] = (*Vphhp)["amie"];
  Fphhp["amie"] += (*Vphpp)["amfe"] * (*Tph)["fi"];

  // II.12
  Fphph["amei"] = (*Vphph)["amei"];
  Fphph["amei"] += (*Vphpp)["amef"] * (*Tph)["fi"];

  // II.13
  Fhpph["ieam"] = (*Vhpph)["ieam"];
  // ERRATUM
  Fhpph["ieam"] += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["an"];

  // II.14
  Fhphp["iema"] = (*Vhphp)["iema"];
  // ERRATUM
  Fhphp["iema"] += (-1.0) * (*Vhhhp)["inme"] * (*Tph)["an"];

  // II.15
  // bug 3
  Fhp["me"] = Whhpp["mnef"] * (*Tph)["fn"];

  // II.1
  Xabie["abie"] = (*Vpphp)["abie"];
  Xabie["abie"] += Fphhh["eimn"] * Xpphh["abnm"];
  Xabie["abie"] += (2.0) * Fphpp["bmef"] * (*Tpphh)["afim"];
  Xabie["abie"] += (-1.0) * Fphpp["bmfe"] * (*Tpphh)["afim"];
  Xabie["abie"] += (-1.0) * Fphpp["bmef"] * (*Tpphh)["afmi"];
  Xabie["abie"] += (-1.0) * Fphpp["amfe"] * (*Tpphh)["bfmi"];
  Xabie["abie"] += (-1.0) * Fphhp["amie"] * (*Tph)["bm"];
  Xabie["abie"] += (-1.0) * Fphph["bmei"] * (*Tph)["am"];
  Xabie["abie"] += (1.0) * (*Vpppp)["abfe"] * (*Tph)["fi"];
  Xabie["abie"] += (-2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfimn"];
  Xabie["abie"] += (1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfnmi"];
  Xabie["abie"] += (1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfinm"];

  // II.2
  Xamij["amij"] = (*Vphhh)["amij"];
  Xamij["amij"] += Fphpp["amef"] * Xpphh["efij"];
  // ERRATUM
  Xamij["amij"] += (2.0) * Fphhh["ejnm"] * (*Tpphh)["aein"];
  // ERRATUM
  Xamij["amij"] += (-1.0) * Fphhh["ejmn"] * (*Tpphh)["aein"];
  // ERRATUM
  Xamij["amij"] += (-1.0) * Fphhh["ejnm"] * (*Tpphh)["eain"];
  // ERRATUM
  Xamij["amij"] += (-1.0) * Fphhh["eimn"] * (*Tpphh)["eajn"];
  Xamij["amij"] += Fhpph["ieam"] * (*Tph)["ej"];
  Xamij["amij"] += Fhphp["jema"] * (*Tph)["ei"];
  Xamij["amij"] += (-1.0) * (*Vhhhh)["ijnm"] * (*Tph)["an"];
  Xamij["amij"] += Fhp["me"] * (*Tpphh)["aeij"];
  Xamij["amij"] += (2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["aefijn"];
  Xamij["amij"] += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["feaijn"];
  Xamij["amij"] += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["afeijn"];

  // II.3
  // NOTE: we antisymmetrize in place Vhphh
  Xim["im"] = (2.0) * (*Vhphh)["iemn"] * (*Tph)["en"];
  // We dont need this term as we have 0 on the lhs
  //      (*Xim)["ii"]       += (*Fepsh)["i"],
  Xim["im"] += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["en"];
  Xim["im"] += Whhpp["mnef"] * Xpphh["efin"];

  // II.4
  // NOTE: we antisymmetrize in place Vphpp
  Xae["ae"] = (2.0) * (*Vphpp)["amef"] * (*Tph)["fm"];
  // We dont need this term as we have 0 on the lhs
  //      (*Xae)["aa"]       += (*Fepsp)["a"],
  Xae["ae"] += (-1.0) * (*Vphpp)["amfe"] * (*Tph)["fm"];
  Xae["ae"] += (-1.0) * Whhpp["mnef"] * Xpphh["afmn"];

  // II.5
  Xjkmn["jkmn"] = (*Vhhhh)["jkmn"];
  // ERRATUM
  Xjkmn["jkmn"] += (*Vhhpp)["mnef"] * Xpphh["efjk"];
  Xjkmn["jkmn"] += (*Vhphh)["jemn"] * (*Tph)["ek"];
  Xjkmn["jkmn"] += (*Vphhh)["ekmn"] * (*Tph)["ej"];

  // II.6
  Xbcef["bcef"] = (*Vpppp)["bcef"];
  // ERRATUM
  Xbcef["bcef"] += (*Vpphh)["efmn"] * Xpphh["bcmn"];
  Xbcef["bcef"] += (-1.0) * (*Vphpp)["bmef"] * (*Tph)["cm"];
  Xbcef["bcef"] += (-1.0) * (*Vhppp)["mcef"] * (*Tph)["bm"];

  // II.7
  Xamei["amei"] = (*Vphph)["amei"];
  Xamei["amei"] += (-1.0) * (*Vhhpp)["mnfe"] * Xpphh["fain"];
  Xamei["amei"] += (-1.0) * (*Vhhph)["nmei"] * (*Tph)["an"];
  Xamei["amei"] += (*Vphpp)["amef"] * (*Tph)["fi"];

  // II.8
  Xamie["amie"] = (*Vphhp)["amie"];
  Xamie["amie"] += Whhpp["mnef"] * (*Tpphh)["afin"];
  Xamie["amie"] += (-1.0) * (*Vhhpp)["mnef"] * Xpphh["fain"];
  Xamie["amie"] += (-1.0) * (*Vhhhp)["nmie"] * (*Tph)["an"];
  Xamie["amie"] += (*Vppph)["aefm"] * (*Tph)["fi"];

  // I.3 until I.10 (numbers are counting terms)
  Xp3h3["abcijk"] = Xjkmn["jkmn"] * (*Tppphhh)["abcimn"];
  Xp3h3["abcijk"] += Xbcef["bcef"] * (*Tppphhh)["aefijk"];
  Xp3h3["abcijk"] += (2.0) * Xamie["amie"] * (*Tppphhh)["ebcmjk"];
  Xp3h3["abcijk"] += (-1.0) * Xamie["amie"] * (*Tppphhh)["becmjk"];
  Xp3h3["abcijk"] += (-1.0) * Xamie["amie"] * (*Tppphhh)["cbemjk"];
  Xp3h3["abcijk"] += (-1.0) * Xamei["amei"] * (*Tppphhh)["ebcmjk"];
  Xp3h3["abcijk"] += (-1.0) * Xamei["bmei"] * (*Tppphhh)["aecmjk"];
  Xp3h3["abcijk"] += (-1.0) * Xamei["cmei"] * (*Tppphhh)["abemjk"];
  Xp3h3["abcijk"] += Xae["ae"] * (*Tppphhh)["ebcijk"];
  Xp3h3["abcijk"] += (-1.0) * Xim["im"] * (*Tppphhh)["abcmjk"];

  // Permute I.3 until I.10 through P(ia/jb, kc)
  (*Rppphhh)["abcijk"] = Xp3h3["abcijk"];
  (*Rppphhh)["abcijk"] += Xp3h3["bacjik"];
  (*Rppphhh)["abcijk"] += Xp3h3["cbakji"];

  // I.1 and I.2
  Xp3h3["abcijk"] = (1.0) * Xabie["abie"] * (*Tpphh)["cekj"];
  Xp3h3["abcijk"] += (-1.0) * Xamij["amij"] * (*Tpphh)["bcmk"];
  (*Rppphhh)["abcijk"] += Xp3h3["abcijk"];
  (*Rppphhh)["abcijk"] += Xp3h3["acbikj"];
  (*Rppphhh)["abcijk"] += Xp3h3["cabkij"];
  (*Rppphhh)["abcijk"] += Xp3h3["cbakji"];
  (*Rppphhh)["abcijk"] += Xp3h3["bcajki"];
  (*Rppphhh)["abcijk"] += Xp3h3["bacjik"];

  // TODO: delete Fij and so on
  return residuum;
}

PTR(FockVector<double>)
CcsdtEnergyFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

PTR(FockVector<sisi4s::complex>)
CcsdtEnergyFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<sisi4s::complex>) &amplitudes) {
  return getResiduumTemplate<sisi4s::complex>(iterationStep, amplitudes);
}
