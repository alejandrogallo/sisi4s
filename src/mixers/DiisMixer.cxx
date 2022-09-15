#include <mixers/DiisMixer.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>
#include <Sisi4s.hpp>
#include <extern/Lapack.hpp>

#include <math/IterativePseudoInverse.hpp>

#include <array>

using namespace CTF;
using namespace sisi4s;

MIXER_REGISTRAR_DEFINITION(DiisMixer);

std::vector<double> inverse(
  std::vector<double> matrix, int N
){
  std::vector<double> column(N,0);
  std::vector<int> ipiv(N);
  std::vector<double> work(N);
  column[0] = -1.0;

  int one(1); int info;
  dsysv_("U", &N, &one, matrix.data(), &N, ipiv.data(), column.data(), &N, work.data(), &N, &info);
  if ( info != 0) {
    LOG(0) << "Problem diagonalization of matrix:\n";
    for (auto i: matrix)  LOG(9) << i << std::endl;
    throw "problem diagonalization\n";
  }
  return column;
}

std::vector<sisi4s::complex> inverse(
  std::vector<sisi4s::complex> matrix, int N
){
  std::vector<sisi4s::complex> column(N);
// TODO

  return column;
}


template <typename F>
DiisMixer<F>::DiisMixer(
  Algorithm *algorithm
):
  Mixer<F>(algorithm), next(nullptr), nextResiduum(nullptr)
{
  int N( algorithm->getRealArgument("maxResidua", 4) );
  LOG(1,"DiisMixer") << "maxResidua=" << N << std::endl;
  EMIT() << YAML::Key << "mixer" << YAML::Value;
  EMIT() << YAML::BeginMap;
  EMIT() << YAML::Key << "type" << YAML::Value << "diis";
  EMIT() << YAML::Key << "max-residua" << YAML::Value << N;
  EMIT() << YAML::EndMap;

  amplitudes.resize(N);
  residua.resize(N);
  nextIndex = 0;
  count = 0;

  // generate initial overlap matrix
  std::array<int,2> lens{{N+1, N+1}};
  std::array<int,2> syms{{NS, NS}};
  B = NEW(CTF::Tensor<F>,
    lens.size(), lens.data(), syms.data(), *Sisi4s::world, "B"
  );
  CTF::Tensor<F> one(false, *B);
  one["ij"] += 1.0;
  std::array<int,2> upperRightBegin{{0,1}};
  std::array<int,2> upperLeftBegin{{1,0}};
  std::array<int,2> upperRightEnd{{1,N+1}};
  std::array<int,2> lowerLeftEnd{{N+1,1}};
  B->slice(
    upperRightBegin.data(), upperRightEnd.data(), 1.0,
    one,
    upperRightBegin.data(), upperRightEnd.data(), -1.0
  );
  B->slice(
    upperLeftBegin.data(), lowerLeftEnd.data(), 1.0,
    one,
    upperLeftBegin.data(), lowerLeftEnd.data(), -1.0
  );
}

template <typename F>
DiisMixer<F>::~DiisMixer() {
}

template <typename F>
void DiisMixer<F>::append(
  const PTR(FockVector<F>) &A, const PTR(FockVector<F>) &R
) {
  // replace amplidue and residuum at nextIndex
  amplitudes[nextIndex] = A;
  residua[nextIndex] = R;

  // generate 1x1 matrix to enter new overlap elements
  const int N(residua.size());
  std::array<int,2> lens{{1, 1}};
  std::array<int,2> syms{{NS, NS}};
  CTF::Tensor<F> one(lens.size(), lens.data(), syms.data(), *Sisi4s::world, "1");
  one["ij"] += 1.0;
  for (int i(0); i < N; ++i) {
    if (residua[i]) {
      std::array<int,2> colBegin{{nextIndex+1,i+1}};
      std::array<int,2> colEnd{{nextIndex+2,i+2}};
      std::array<int,2> rowBegin{{i+1,nextIndex+1}};
      std::array<int,2> rowEnd{{i+2,nextIndex+2}};
      std::array<int,2> oneBegin{{0,0}};
      std::array<int,2> oneEnd{{1,1}};
      F overlap( 2.0*std::real(residua[i]->dot(*R)) );
      B->slice(
        colBegin.data(), colEnd.data(), 0.0,
        one,
         oneBegin.data(), oneEnd.data(), overlap
      );
      if (i == nextIndex) continue;
      B->slice(
        rowBegin.data(), rowEnd.data(), 0.0,
        one,
        oneBegin.data(), oneEnd.data(), overlap
      );
    }
  }
  if (count < N) ++count;

//  B->print();
  // now, pseudo-invert upper left corner of B and read out its first column
  int dim(count+1);
  std::array<int,2> upperLeftBegin{{0, 0}};
  std::array<int,2> lowerRightEnd{{count+1, count+1}};
  std::array<int,2> firstColEnd{{count+1,1}};
  std::vector<F> column(count+1);
  std::vector<F> matrix(dim*dim);
  B->slice(upperLeftBegin.data(), lowerRightEnd.data()).read_all(matrix.data());
  column = inverse(matrix, dim);

  LOG(1, "DiisMixer") << "lambda" << "=" << column[0] << std::endl;
  EMIT() << YAML::Key << "mixing";
  EMIT() << YAML::Value << YAML::BeginMap;
  // FIXME: imaginary part dismissed
  EMIT() << YAML::Key << "lambda" << YAML::Value << std::real(column[0]);

  next = NEW(FockVector<F>, *A);
  *next *= F(0);
  nextResiduum = NEW(FockVector<F>, *R);
  *nextResiduum *= F(0);
  EMIT() << YAML::Key << "weights";
  EMIT() << YAML::Value << YAML::Flow << YAML::BeginSeq;
  for (int j(0); j < count; ++j) {
    int i( (nextIndex+N-j) % N );
    LOG(1, "DiisMixer") << "w^(-" << (j+1) << ")=" << column[i+1] << std::endl;
    // FIXME: imaginary part dismissed
    EMIT() << std::real(column[i+1]);
    *next += column[i+1] * *amplitudes[i];
    *nextResiduum += column[i+1] * *residua[i];
  }
  EMIT() << YAML::EndSeq;
  EMIT() << YAML::Comment("w(-1), ... , w(-max-residua)");
  EMIT() << YAML::EndMap;

  nextIndex = (nextIndex+1) % N;
}

template <typename F>
PTR(const FockVector<F>) DiisMixer<F>::get() {
  return next;
}

template <typename F>
PTR(const FockVector<F>) DiisMixer<F>::getResiduum() {
  return nextResiduum;
}


// instantiate
template class sisi4s::DiisMixer<sisi4s::Float64>;
template class sisi4s::DiisMixer<sisi4s::Complex64>;

