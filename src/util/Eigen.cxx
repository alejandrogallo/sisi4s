#include <util/Eigen.hpp>
#include <Sisi4s.hpp>

namespace sisi4s {

MatrixColumnMajor toEigenMatrix(Tensor<double> &ctf) {
  MatrixColumnMajor result(ctf.lens[0], ctf.lens[1]);

  const size_t size(ctf.lens[0] * ctf.lens[1]);
  std::vector<int64_t> indices(size);
  std::vector<double> values(size);

  // make indices a range from 0 to indices.size()-1
  std::iota(indices.begin(), indices.end(), 0);

  ctf.read(indices.size(), indices.data(), values.data());

  for (size_t i(0); i < values.size(); i++) { result(i) = values[i]; }

  return result;
}

Tensor<double> toCtfMatrix(const MatrixColumnMajor &m) {
  int syms[] = {NS, NS}, lens[] = {(int)m.rows(), (int)m.cols()};
  const int64_t size(m.rows() * m.cols()), r(Sisi4s::world->rank),
      np(Sisi4s::world->np), chunks(size / np), start(r * chunks),
      end((np - 1) == r ? size : (r + 1) * chunks);
  ;
  std::vector<int64_t> indices(end - start);
  Tensor<double> t(2, lens, syms, *Sisi4s::world);

  //  make indices a range from 0 to indices.size()-1
  std::iota(indices.begin(), indices.end(), start);
  t.write(indices.size(), indices.data(), &m(0) + start);

  return t;
}

} // namespace sisi4s
