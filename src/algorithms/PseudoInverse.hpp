/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PSEUDO_INVERSE_DEFINED
#define PSEUDO_INVERSE_DEFINED

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
class PseudoInverse : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(PseudoInverse);
  PseudoInverse(std::vector<Argument> const &argumentList);
  virtual ~PseudoInverse();
  virtual void run();
};
} // namespace sisi4s

#endif
