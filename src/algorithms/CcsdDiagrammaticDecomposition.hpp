/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_DIAGRAMMATIC_DECOMPOSITION
#define CCSD_DIAGRAMMATIC_DECOMPOSITION

#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <vector>
#include <string>

namespace sisi4s {
  class CcsdDiagrammaticDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdDiagrammaticDecomposition);
    CcsdDiagrammaticDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcsdDiagrammaticDecomposition();
    virtual void run();
  protected:
    void evaluateEnergy(std::string diagramType,
                        Tensor<double> &deltaabij,
                        Tensor<double> &deltaai,
                        Tensor<double> &Rabij,
                        Tensor<double> &Rai);
   void sliceIntoResiduum(Tensor<double> &Rxyij,
                          int a,
                          int b,
                          Tensor<double> &Rabij);
  };
}
#endif
