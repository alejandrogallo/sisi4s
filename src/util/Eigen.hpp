#ifndef EIGEN_HPP_
#define EIGEN_HPP_

#include <Eigen/Dense>
#include <util/Tensor.hpp>

namespace sisi4s {

using MatrixColumnMajor =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

MatrixColumnMajor toEigenMatrix(Tensor<double> &ctf);
Tensor<double> toCtfMatrix(const MatrixColumnMajor &m);

} // namespace sisi4s

#endif
