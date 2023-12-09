#include <string>
#include <vector>
#include <algorithm>
#include <util/Libint.hpp>
#include <algorithms/OneBodyFromGaussian.hpp>
#include <util/Tensor.hpp>
#include <Sisi4s.hpp>
#include <util/Log.hpp>
#include <util/Integrals.hpp>
#include <iostream>
#include <Eigen/Eigen>
#include <numeric>
#include <set>
#include <map>
#include <util/Emitter.hpp>
#define LOGGER(_l) LOG(_l, "OneBodyFromGaussian")

using namespace sisi4s;

using RowMajor =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

IMPLEMENT_EMPTY_DRYRUN(OneBodyFromGaussian) {}

struct Sinfo {
  const size_t size, begin;
  inline Sinfo(const libint2::BasisSet &shells, const size_t i)
      : size(shells[i].size())
      , begin(shells.shell2bf()[i]) {}
};

static RowMajor getOneBodyIntegrals(const libint2::BasisSet &shells,
                                    const libint2::Operator op,
                                    const std::vector<libint2::Atom> &atoms) {

  // Get number of basis set functions
  const size_t Np = shells.nbf();
  RowMajor result(RowMajor::Zero(Np, Np));

  // construct the overlap integrals engine
  libint2::Engine engine(op, shells.max_nprim(), shells.max_l(), 0);

  if (op == libint2::Operator::nuclear) {
    engine.set_params(libint2::make_point_charges(atoms));
  }

  const auto &resultBuffer = engine.results();

  // loop over unique shell pairs p>=q
  for (size_t p = 0; p < shells.size(); ++p) {
    for (size_t q = 0; q <= p; ++q) {

      const Sinfo P(shells, p), Q(shells, q);

      engine.compute(shells[p], shells[q]);

      // copy the result buffer of size (P.size*Q.size) to the
      // bufferMatrix, needed to do a block assignment to the Eigen matrix
      Eigen::Map<const RowMajor> bufferMatrix(resultBuffer[0], P.size, Q.size);
      result.block(P.begin, Q.begin, P.size, Q.size) = bufferMatrix;

      if (p != q) {
        // transpose and copy to the non-explicitly calculated block
        result.block(Q.begin, P.begin, Q.size, P.size) =
            bufferMatrix.transpose();
      }

    } // q
  }   // p

  return result;
}

DEFSPEC(OneBodyFromGaussian,
        SPEC_IN({"basisSet",
                 SPEC_VALUE("The name of the basis set", std::string)},
                {"kernel",
                 SPEC_ONE_OF("The type of 1-body integral to evaluate",
                             std::string,
                             "nuclear",
                             "overlap",
                             "kinetic")},
                {"atoms",
                 SPEC_VARIN(
                     "Vector of libint atoms specifying a molecular structure",
                     std::vector<libint2::Atom> *)}),
        SPEC_OUT({"Out", SPEC_VAROUT("TODO: DOC", Tensor<double> *)}));

IMPLEMENT_ALGORITHM(OneBodyFromGaussian) {

  libint2::initialize();

  const std::string basisSet(in.get<std::string>("basisSet")),
      kernel(in.get<std::string>("kernel"));

  auto const &atoms = *in.get<std::vector<libint2::Atom> *>("atoms");
  const libint2::BasisSet shells(basisSet, atoms);
  const int Np(shells.nbf());
  const libint2::Operator op([&](void) {
    if (kernel == "nuclear") return libint2::Operator::nuclear;
    else if (kernel == "overlap") return libint2::Operator::overlap;
    else if (kernel == "kinetic") return libint2::Operator::kinetic;
    else throw EXCEPTION("Kernel not recognised");
  }());

  LOGGER(1) << "kernel: " << kernel << std::endl;

  EMIT() << YAML::Key << "Np" << YAML::Value << Np << YAML::Key << "basis-set"
         << YAML::Value << basisSet << YAML::Key << "kernel" << YAML::Value
         << kernel;

  auto result(getOneBodyIntegrals(shells, op, atoms));
  double *data(&result(0));
  LOGGER(1) << "computation done" << std::endl;
  libint2::finalize();

  const int rank_m = int(Sisi4s::world->rank == 0); // rank mask
  const std::vector<int> lens(2, Np);
  const std::vector<int> syms(2, NS);
  auto Opq(new Tensor<double>(2, lens.data(), syms.data(), *Sisi4s::world));

  {
    std::vector<int64_t> indices(rank_m * Np * Np);
    std::iota(indices.begin(), indices.end(), 0);
    LOGGER(1) << "writing into ctf tensor" << std::endl;
    Opq->write(indices.size(), indices.data(), data);
  }

  out.set<Tensor<double> *>("Out", Opq);
}
