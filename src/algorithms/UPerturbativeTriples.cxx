#include <algorithms/UPerturbativeTriples.hpp>
#include <math/MathFunctions.hpp>
#include <DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Sisi4s.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

ALGORITHM_REGISTRAR_DEFINITION(UPerturbativeTriples);

UPerturbativeTriples::UPerturbativeTriples(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

UPerturbativeTriples::~UPerturbativeTriples() {
}

void UPerturbativeTriples::run() {
  Tensor<double>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<double>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<double> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<double> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<double> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<double> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<double>   *Tai(getTensorArgument("CcsdSinglesAmplitudes"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  int vvvooo[] = { Nv, Nv, Nv, No, No, No };
  int   syms[] = { NS, NS, NS, NS, NS, NS };
  Tensor<double> SVabcijk(6, vvvooo, syms, *Vabij->wrld, "SVabcijk");
  // TODO: erase 0.5
  SVabcijk["abcijk"]  = 0.5 * (*Tai)["ai"] * (*Vabij)["bcjk"];

  Tensor<double> DVabcijk(6, vvvooo, syms, *Vabij->wrld, "DVabcijk");
  DVabcijk["abcijk"]  =          (*Tabij)["adij"] * (*Vabci)["bcdk"];
  DVabcijk["abcijk"] += (-1.0) * (*Tabij)["abil"] * (*Vijka)["jklc"];

  Tensor<double> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");

  /* AUTOMATIC
               Tabij["adij"]   Vabci["bcdk"]
  ( + 1.0  ) * Tabij["dgjk"] * Vabci["efgi"];  ⇒ defjki
  ( - 1.0  ) * Tabij["egjk"] * Vabci["dfgi"];  ⇒ edfjki
  ( - 1.0  ) * Tabij["fgjk"] * Vabci["edgi"];  ⇒ fedjki
  ( + 1.0  ) * Tabij["dgki"] * Vabci["efgj"];  ⇒ defkij
  ( - 1.0  ) * Tabij["egki"] * Vabci["dfgj"];  ⇒ edfkij
  ( - 1.0  ) * Tabij["fgki"] * Vabci["edgj"];  ⇒ fedkij
  ( + 1.0  ) * Tabij["dgij"] * Vabci["efgk"];  ⇒ defijk
  ( - 1.0  ) * Tabij["egij"] * Vabci["dfgk"];  ⇒ edfijk
  ( - 1.0  ) * Tabij["fgij"] * Vabci["edgk"];  ⇒ fedijk

  AUTOMATIC HOLES

             - Tabij["abli"] * Vijka["jklc"]; (ours)
  ( + 1.0  ) * Tabij["deok"] * Vijka["ijof"];  ⇒ defkij
  ( - 1.0  ) * Tabij["dfok"] * Vijka["ijoe"];  ⇒ dfekij
  ( - 1.0  ) * Tabij["feok"] * Vijka["ijod"];  ⇒ fedkij
  ( - 1.0  ) * Tabij["deoj"] * Vijka["ikof"];  ⇒ defjik
  ( + 1.0  ) * Tabij["dfoj"] * Vijka["ikoe"];  ⇒ dfejik
  ( + 1.0  ) * Tabij["feoj"] * Vijka["ikod"];  ⇒ fedjik
  ( - 1.0  ) * Tabij["deoi"] * Vijka["kjof"];  ⇒ defikj
  ( + 1.0  ) * Tabij["dfoi"] * Vijka["kjoe"];  ⇒ dfeikj
  ( + 1.0  ) * Tabij["feoi"] * Vijka["kjod"];  ⇒ fedikj

  AUTOMATIC singles

  ( + 1.0  ) * Tai["dk"] * Vabij["efij"];defkij
  ( - 1.0  ) * Tai["ek"] * Vabij["dfij"];edfkij
  ( - 1.0  ) * Tai["fk"] * Vabij["edij"];fedkij
  ( - 1.0  ) * Tai["dj"] * Vabij["efik"];defjik
  ( + 1.0  ) * Tai["ej"] * Vabij["dfik"];edfjik
  ( + 1.0  ) * Tai["fj"] * Vabij["edik"];fedjik
  ( - 1.0  ) * Tai["di"] * Vabij["efkj"];defikj
  ( + 1.0  ) * Tai["ei"] * Vabij["dfkj"];edfikj
  ( + 1.0  ) * Tai["fi"] * Vabij["edkj"];fedikj

  */

  /*
   * Explanation
   * -----------
   * We have to anisymmetrize the diagramas because we are doing
   * unrestricted and hier everything is a matrix element.
   *
   * We can't put the VIJKA diagram with the VABCI diagram
   * because they don't antisymmetrize in the same way, think about it,
   * or just draw them to see.
   *
   * So first antisymmetrize the VABCI-based diagram.
   * Then we antisymmetrize the VIJKA-based diagram
   *
   * We store the sum of these antisymmetrized diagram into DVabcijk
   * and we still add the singles contribution to Tabcijk.
   *
   * We divide Tabcijk by the denominator and calculate the energy.
   *
   */

  Tabcijk["defjki"]  = 0.0;
  SVabcijk["abcijk"] = 0.0;

  // VABCI PART   -------------------------------------------------------------
  DVabcijk["abcijk"] = (*Tabij)["adij"] * (*Vabci)["bcdk"];
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defjki"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfjki"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedjki"];
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defkij"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfkij"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedkij"];
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defijk"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfijk"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedijk"];

  // SAVE antisymmetrized VABCI PART FOR LATER
  //SVabcijk["abcijk"] += Tabcijk["abcijk"];

  // VIJKA PART   -------------------------------------------------------------
  DVabcijk["defkij"] = (*Tabij)["deok"] * (*Vijka)["ijof"];
  //--
  Tabcijk["defjki"] += ( + 1.0  ) * DVabcijk["defkij"];
  Tabcijk["defjki"] += ( - 1.0  ) * DVabcijk["dfekij"];
  Tabcijk["defjki"] += ( - 1.0  ) * DVabcijk["fedkij"];
  Tabcijk["defjki"] += ( - 1.0  ) * DVabcijk["defjik"];
  Tabcijk["defjki"] += ( + 1.0  ) * DVabcijk["dfejik"];
  Tabcijk["defjki"] += ( + 1.0  ) * DVabcijk["fedjik"];
  Tabcijk["defjki"] += ( - 1.0  ) * DVabcijk["defikj"];
  Tabcijk["defjki"] += ( + 1.0  ) * DVabcijk["dfeikj"];
  Tabcijk["defjki"] += ( + 1.0  ) * DVabcijk["fedikj"];

  // Add antisymmetrized VIJKA part
  //SVabcijk["abcijk"] += Tabcijk["abcijk"];

  // Save antisymmetrized in DVabcijk
  DVabcijk["abcijk"] = Tabcijk["abcijk"];

  // Singles part -------------------------------------------------------------
  SVabcijk["defkij"] = ( + 1.0  ) * (*Tai)["dk"] * (*Vabij)["efij"];
  // --
  Tabcijk["defkij"] += ( + 1.0  ) * SVabcijk["defkij"];
  Tabcijk["defkij"] += ( - 1.0  ) * SVabcijk["edfkij"];
  Tabcijk["defkij"] += ( - 1.0  ) * SVabcijk["fedkij"];
  Tabcijk["defkij"] += ( - 1.0  ) * SVabcijk["defjik"];
  Tabcijk["defkij"] += ( + 1.0  ) * SVabcijk["edfjik"];
  Tabcijk["defkij"] += ( + 1.0  ) * SVabcijk["fedjik"];
  Tabcijk["defkij"] += ( - 1.0  ) * SVabcijk["defikj"];
  Tabcijk["defkij"] += ( + 1.0  ) * SVabcijk["edfikj"];
  Tabcijk["defkij"] += ( + 1.0  ) * SVabcijk["fedikj"];


  // HIRATA Implementation
  // ---------------------
  //
  // Let us find index equivalences between hirata's and ours
  //
  //           Tabij["gdjk"]    Vabic["efig"]   (original)
  //           Tabij["dgjk"]    Vabci["efgi"]   (changed for our code)
  //           Tabij["adij"]    Vabci["bcdk"]   (ours)
  //           => defjki
  // and
  //
  //           Tabij["deok"] * Vijka["ijof"];
  //         - Tabij["abli"] * Vijka["jklc"];
  //           => defjki
  // and
  //          Tai["dk"] * Vabij["efij"];
  //          Tai["ai"] * Vabij["bcjk"];
  //          => defjki

#ifdef OLD_HIRATA
  Tabcijk["defjki"] = 0.0;
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defjki"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfjki"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedjki"];
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defkij"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfkij"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedkij"];
  //--
  Tabcijk["defjki"] += (+ 1.0) * DVabcijk["defijk"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["edfijk"];
  Tabcijk["defjki"] += (- 1.0) * DVabcijk["fedijk"];

  //
  DVabcijk["defjki"] =  Tabcijk["defjki"];

  /*
  //--
    ( + 1.0  ) * Tabij["gdjk"] * Vabic["efig"];
    ( - 1.0  ) * Tabij["gejk"] * Vabic["dfig"];
    ( - 1.0  ) * Tabij["gfjk"] * Vabic["edig"];
  //--
    ( + 1.0  ) * Tabij["gdki"] * Vabic["efjg"];
    ( - 1.0  ) * Tabij["geki"] * Vabic["dfjg"];
    ( - 1.0  ) * Tabij["gfki"] * Vabic["edjg"];
  //--
    ( + 1.0  ) * Tabij["gdij"] * Vabic["efkg"];
    ( - 1.0  ) * Tabij["geij"] * Vabic["dfkg"];
    ( - 1.0  ) * Tabij["gfij"] * Vabic["edkg"];
    */

  Tabcijk["defjki"] += (+ 1.0) * SVabcijk["defjki"];
  Tabcijk["defjki"] += (- 1.0) * SVabcijk["edfjki"];
  Tabcijk["defjki"] += (- 1.0) * SVabcijk["fedjki"];
  Tabcijk["defjki"] += (- 1.0) * SVabcijk["defkji"];
  Tabcijk["defjki"] += (+ 1.0) * SVabcijk["edfkji"];
  Tabcijk["defjki"] += (+ 1.0) * SVabcijk["fedkji"];
  Tabcijk["defjki"] += (- 1.0) * SVabcijk["defjik"];
  Tabcijk["defjki"] += (+ 1.0) * SVabcijk["edfjik"];
  Tabcijk["defjki"] += (+ 1.0) * SVabcijk["fedjik"];
#endif

  //( + 1.0  ) * Tai["dk"] * Vabij["efij"];
  //( - 1.0  ) * Tai["ek"] * Vabij["dfij"];
  //( - 1.0  ) * Tai["fk"] * Vabij["edij"];
  //( - 1.0  ) * Tai["dj"] * Vabij["efik"];
  //( + 1.0  ) * Tai["ej"] * Vabij["dfik"];
  //( + 1.0  ) * Tai["fj"] * Vabij["edik"];
  //( - 1.0  ) * Tai["di"] * Vabij["efkj"];
  //( + 1.0  ) * Tai["ei"] * Vabij["dfkj"];
  //( + 1.0  ) * Tai["fi"] * Vabij["edkj"];



  /**********************
  // Gauss implementaiton
  // permute cylicly abc
  //    permute cyclicly ijk
  Tabcijk["abcijk"]  = (+1.0) * DVabcijk["abcijk"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["abcjki"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["abckij"];
  //    permute cyclicly ijk
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["bcaijk"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["bcajki"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["bcakij"];
  //    permute cyclicly ijk
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["cabijk"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["cabjki"];
  Tabcijk["abcijk"] += (+1.0) * DVabcijk["cabkij"];
  */

  /**********************
  // Gauss implementaiton
  // permute cylicly abc
  //    permute cyclicly ijk
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["abcijk"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["abcjki"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["abckij"];
  //    permute cyclicly ijk
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["bcaijk"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["bcajki"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["bcakij"];
  //    permute cyclicly ijk
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["cabijk"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["cabjki"];
  Tabcijk["abcijk"] += (+1.0) * SVabcijk["cabkij"];
  */

  /***********************
  // Theo's implementation
  //  Tabcijk["abcijk"] += (+1.0) * DVabcijk["acbijk"];
  //  Tabcijk["abcijk"] += (+1.0) * DVabcijk["bacijk"];
  //  Tabcijk["abcijk"] += (+1.0) * DVabcijk["bcaijk"];
  //  Tabcijk["abcijk"] += (+1.0) * DVabcijk["cabijk"];
  //  Tabcijk["abcijk"] += (+1.0) * DVabcijk["cbaijk"];
  */

  /***********************
  // Theo's implementation
  //Tabcijk["abcijk"] += (+1.0) * SVabcijk["acbijk"];
  //Tabcijk["abcijk"] += (+1.0) * SVabcijk["bacijk"];
  //Tabcijk["abcijk"] += (+1.0) * SVabcijk["bcaijk"];
  //Tabcijk["abcijk"] += (+1.0) * SVabcijk["cabijk"];
  //Tabcijk["abcijk"] += (+1.0) * SVabcijk["cbaijk"];
  */

  SVabcijk["abcijk"]  =          (*epsi)["i"];
  SVabcijk["abcijk"] +=          (*epsi)["j"];
  SVabcijk["abcijk"] +=          (*epsi)["k"];
  SVabcijk["abcijk"] += (-1.0) * (*epsa)["a"];
  SVabcijk["abcijk"] += (-1.0) * (*epsa)["b"];
  SVabcijk["abcijk"] += (-1.0) * (*epsa)["c"];
  Bivar_Function<> fDivide(&divide<double>);
  Tabcijk.contract(
    1.0, Tabcijk,"abcijk", SVabcijk,"abcijk", 0.0,"abcijk", fDivide
  );

  Scalar<> energy(*Sisi4s::world);
  energy[""]  = (1.0 / 36.0) * DVabcijk["abcijk"] * Tabcijk["abcijk"];
//  energy[""] += DVabcijk["bacjik"] * Tabcijk["abcijk"];
//  energy[""] += DVabcijk["acbikj"] * Tabcijk["abcijk"];
//  energy[""] += DVabcijk["cbakji"] * Tabcijk["abcijk"];
//  energy[""] += DVabcijk["cabkij"] * Tabcijk["abcijk"];
//  energy[""] += DVabcijk["bcajki"] * Tabcijk["abcijk"];

  double eTriples(energy.get_val());
  //double eCcsd(getRealArgument("CcsdEnergy"));
  //double e(eCcsd + eTriples);
  //LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  //LOG(1, "PerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  
  LOG(0, "PerturbativeTriples") << "triples=" << eTriples << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", energy);
}

void UPerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PPPHCoulombIntegrals");

  DryTensor<> *Tai(
    getTensorArgument<double, DryTensor<double>>("CcsdSinglesAmplitudes")
  );
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("CcsdDoublesAmplitudes")
  );

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int vvvooo[] = { Nv, Nv , Nv , No , No , No };
  int   syms[] = { NS, NS,  NS , NS , NS , NS };
  DryTensor<> Tabcijk(6, vvvooo, syms, SOURCE_LOCATION);

  {
    DryTensor<> Zabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  }

  DryTensor<> Zai(*Tai, SOURCE_LOCATION);
  DryTensor<> Zabij(*Tabij, SOURCE_LOCATION);

  DryScalar<> energy();
}


