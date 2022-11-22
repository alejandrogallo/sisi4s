/* Copyright 2020. Erwin Schroedinger and Egyl Hylleraas */
#ifndef __MEAN_CORRELATION_HOLE_DEPTH_BLAH_DEFINED

#  include <algorithms/Algorithm.hpp>

namespace sisi4s {
struct MeanCorrelationHoleDepth : public Algorithm {
  ALGORITHM_REGISTRAR_DECLARATION(MeanCorrelationHoleDepth);

  MeanCorrelationHoleDepth(std::vector<Argument> const &argumentList)
      : Algorithm(argumentList) {}
  ~MeanCorrelationHoleDepth(){};

  virtual void run();
};
} // namespace sisi4s

#endif
