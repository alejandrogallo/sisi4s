#ifndef _HARTREE_FOCK_FROM_COULOMB_INTEGRALS_ALGO_DEFINED
#define _HARTREE_FOCK_FROM_COULOMB_INTEGRALS_ALGO_DEFINED

#include <algorithms/Algorithm.hpp>
#include <Eigen/Dense>
#include <util/Tensor.hpp>

namespace sisi4s {

using MatrixColumnMajor =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

class HartreeFockFromCoulombIntegrals : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(HartreeFockFromCoulombIntegrals);
  HartreeFockFromCoulombIntegrals(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~HartreeFockFromCoulombIntegrals() {}

  virtual void run();
};

MatrixColumnMajor toEigenMatrix(Tensor<double> &ctf);

} // namespace sisi4s

#endif
