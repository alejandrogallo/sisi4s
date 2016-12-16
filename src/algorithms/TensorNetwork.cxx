#include <algorithms/TensorNetwork.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <Cc4s.hpp>

#include <vector>
#include <memory>
using std::shared_ptr;

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(TensorNetwork);

TensorNetwork::TensorNetwork(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorNetwork::~TensorNetwork() {
}

/**
 * \brief Testing environement
 */
void TensorNetwork::run() {
}


void TensorNetwork::dryRun() {
  int No(10);
  int Nv(90);
  int Np(No+Nv);
  int NF(200);
  int NR(300);
  auto machineTensorFactory(
//    CtfMachineTensorFactory<>::create(Cc4s::world)
    DryMachineTensorFactory<>::create()
  );
  auto tcc(Tcc<>::create(machineTensorFactory));

/*
  shared_ptr<Tensor<complex>> Tc(
    tcc.createTensor<complex>(std::vector<int>({100,100,10,10}), "Tc")
  );
*/
  auto T(
    tcc->createTensor(std::vector<int>({100,100,10,10}), "T")
  );
  auto Pi(
    tcc->createTensor(std::vector<int>({300,100}), "Pi")
  );
  auto PiT(
    tcc->createTensor(std::vector<int>({300,100}), "PiT")
  );
  auto Lambda(
    tcc->createTensor(std::vector<int>({300,200}), "Lambda")
  );
  auto LambdaT(
    tcc->createTensor(std::vector<int>({300,200}), "LambdaT")
  );

//  CompoundDryTensorExpression<> Gamma("Fac") = PiT["Ra"] * Pi["Rc"] * Lambda["RG"]

  auto ladderOperation = tcc->compile(
    (
      (*T)["abij"] -= -1/4. *(*T)["abji"],
      (*T)["abij"] <<=
        2 * (*T)["cdij"] * (*Pi)["Rd"] * (*PiT)["Rb"] *
        (*Pi)["Sc"] * (*PiT)["Sa"] * (-3) * (*LambdaT)["SF"] * (*Lambda)["RF"]
    )
  );
  ladderOperation->execute();

// this contraction already requires heuristics
/*
  shared_ptr<Tensor<>> Pia(
    tcc.createTensor<>(std::vector<int>({NR,Nv}), "Pia")
  );
  shared_ptr<Tensor<>> Pii(
    tcc.createTensor<>(std::vector<int>({NR,No}), "Pii")
  );
  int Nn(7);
  shared_ptr<Tensor<>> w(
    tcc.createTensor<>(std::vector<int>({Nn}), "w")
  );
  shared_ptr<Tensor<>> H(
    tcc.createTensor<>(std::vector<int>({No,Nn}), "H")
  );
  shared_ptr<Tensor<>> P(
    tcc.createTensor<>(std::vector<int>({Nv,Nn}), "P")
  );
  shared_ptr<Tensor<>> e(
    tcc.createTensor<>(std::vector<int>(), "e")
  );

  shared_ptr<Operation<>> imaginaryTimeMp2Operation = compile(
    (*e)[""] <<=
      (*Pii)["Ri"]  * (*Pia)["Ra"] *
        (*LambdaT)["RF"] * (*Lambda)["SF"] *
      (*Pii)["Sj"] * (*Pia)["Sb"] *
        (*w)["n"] * (*P)["an"] * (*H)["in"] * (*P)["bn"] * (*H)["jn"] *
      (*Pii)["Ti"]  * (*Pia)["Ta"] *
        (*LambdaT)["TH"] * (*Lambda)["UH"] *
      (*Pii)["Uj"] * (*Pia)["Ub"]
  );
  imaginaryTimeMp2Operation->execute();
*/
}


