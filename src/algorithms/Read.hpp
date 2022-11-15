#ifndef READ_HPP_
#define READ_HPP_
#include <algorithms/Algorithm.hpp>

namespace sisi4s {

  namespace cc4s {

    enum AxisType {
      AuxiliaryField,
      State
    };

    struct Dimension {
      size_t length;
      AxisType type;
    };
    using Dimensions = std::vector<Dimension>;

    enum ElementFileType {
      TextFile,
      IeeeBinaryFile
    };

    enum ScalarType {
      Real64,
      Complex64
    };

    template <ScalarType t>
    struct ScalarTypeTraits;

    template <>
    struct ScalarTypeTraits<ScalarType::Real64> {
      using type = double;
    };

    template <>
    struct ScalarTypeTraits<ScalarType::Complex64> {
      using type = sisi4s::complex;
    };

    enum ReadableType {
      Tensor,
    };

    struct ReadHeader {

      enum Version {
        ONE = 100
      };

      Version version;
      ReadableType type;
      ScalarType scalarType;
      Dimensions dimensions;
      ElementFileType elementsType;
      double unit;
    };

  }  // cc4s

  DEFINE_ALGORITHM_HEADER(Read,);
  DEFINE_ALGORITHM_HEADER(Write,);

}  // namespace sisi4s

#endif
