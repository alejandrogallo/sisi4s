#ifndef _ONE_BODY_FROM_GAUSSIAN_DEFINED
#define _ONE_BODY_FROM_GAUSSIAN_DEFINED

#include <vector>

#include <algorithms/Algorithm.hpp>
#include <math/Complex.hpp>
#include <util/Libint.hpp>
#include <Eigen/Eigen>

namespace sisi4s {

class OneBodyFromGaussian : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(OneBodyFromGaussian);
  OneBodyFromGaussian(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~OneBodyFromGaussian() {}
  virtual void run();
};

} // namespace sisi4s

#endif
