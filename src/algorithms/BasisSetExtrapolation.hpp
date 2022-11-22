/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED
#define BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED

#include <vector>
#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {
class BasisSetExtrapolation : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(BasisSetExtrapolation);
  BasisSetExtrapolation(std::vector<Argument> const &argumentList);
  virtual ~BasisSetExtrapolation();
  virtual void run();

protected:
  void evaluateQGG(int orbitalPairStart, int orbtialPairEnd, int slice);
  void fitF12(int type, real minG, real maxG);
  void calculateNewSF(int type,
                      real gamma,
                      Tensor<double> *coulombKernel,
                      Tensor<double> *newSF,
                      Tensor<double> *resNewSF);
  void invertQGG();
};
} // namespace sisi4s

#endif
