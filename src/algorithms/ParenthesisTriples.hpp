#ifndef PARENTHESIS_TRIPLES
#define PARENTHESIS_TRIPLES

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
/**
 * \brief Caclulates perturbative triples correction, third attemp
 */
struct IJKPointer {
  double *i, *j, *k;
};

template <typename F>
struct IrmlerTensor {
  std::vector<int64_t> lens;
  int64_t order;
  std::vector<F> data;
  // A move constructor should be created automatically by the compiler
  IrmlerTensor(std::vector<int64_t> lens_, int64_t order_, std::vector<F> data_)
      : lens(lens_)
      , order(order_)
      , data(data_) {}
};

template <typename F>
struct VectorTensor {
  static const int64_t NO_REF = -1;
  std::array<IrmlerTensor<F> *, 3> tensors;
  std::array<int64_t, 3> refs;
  IrmlerTensor<F> *withReference(int64_t r) {
    auto *p = std::find(refs.begin(), refs.end(), r);
    //      if (p == refs.size()) {/* panic */}
    const int64_t idx(std::distance(refs.begin(), p));
    return tensors[idx];
  }
  VectorTensor()
      : tensors({nullptr, nullptr, nullptr})
      , refs({NO_REF, NO_REF, NO_REF}) {}
};

DEFINE_ALGORITHM_HEADER(

    ParenthesisTriples,

    int64_t No;
    int64_t Nv;
    int64_t NG;
    int64_t NvCube;
    double *Tabc;
    double *Zabc;
    double *scratch;
    double *Tabij, *Tai;
    double *epsi, *epsa;
    double *Vabij, *Vhphh, *Vppph = nullptr;
    double *Vpppijk = nullptr;
    double *realGab = nullptr;
    double *imagGab = nullptr;
    double *realGai = nullptr;
    double *imagGai = nullptr;
    bool particleDiagram = true, holeDiagram = true;
    double getEnergy(const double epsijk,
                     const double *Tabc_,
                     const double *Zabc_);
    void doublesContribution(const std::array<int64_t, 3> &ijk,
                             IJKPointer integralContainer,
                             double *scratchV,
                             double *scratchO,
                             double *output);
    void singlesContribution(const std::array<int64_t, 3> &ijk,
                             double *scratchO,
                             double *output);
    void getVpppijkFromVertex(const std::array<int64_t, 3> &ijk,
                              double *scratchO,
                              double *output);

    IJKPointer getVpppijkOnTheFly(const std::array<int64_t, 3> &ijk,
                                  VectorTensor<double> &vtensor););
} // namespace sisi4s

#endif
