#ifndef UEG_VERTEX_GENERATOR_DEFINED
#define UEG_VERTEX_GENERATOR_DEFINED

#include <array>

#include <algorithms/Algorithm.hpp>

namespace sisi4s {
using ivec = std::array<int, 3>;
using dvec = std::array<double, 4>;

class UegVertexGenerator : public Algorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UegVertexGenerator);
  ~UegVertexGenerator();
  UegVertexGenerator(std::vector<Argument> const &argumentList);
  void run();

  template <typename F>
  void run();

protected:
  double evalMadelung(double volume);
  double Vijji(const dvec a, const dvec b, const double v);

  bool halfGrid;
  size_t No, Nv, NF;
  double rs, madelung;
};
} // namespace sisi4s

#endif
