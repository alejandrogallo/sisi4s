#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <util/Scanner.hpp>
#include <util/Tensor.hpp>

namespace sisi4s {

namespace cc4s {

enum AxisType { AuxiliaryField, State };

struct Dimension {
  size_t length;
  AxisType type;
};
using Dimensions = std::vector<Dimension>;

enum ElementFileType { TextFile, IeeeBinaryFile };

enum ScalarType { Real64, Complex64 };

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

  enum Version { ONE = 100 };

  Version version;
  ReadableType type;
  ScalarType scalarType;
  Dimensions dimensions;
  ElementFileType elementsType;
  double unit;
};

} // namespace cc4s

class TensorIo {
public:
  template <typename F = real, typename T = Tensor<F>>
  static T *readBinary(std::string const &fileName);
  template <typename F = real, typename T = Tensor<F>>
  static T *readText(std::string const &fileName,
                     std::string const &delimiter = " ",
                     int64_t const bufferSize = 1024 * 1024 * 1024);

  template <typename F = real, typename T = Tensor<F>>
  static void writeBinary(std::string const &fileName, T &A);
  template <typename F = real, typename T = Tensor<F>>
  static void writeText(std::string const &fileName,
                        T &A,
                        std::string const &rowIndexOrder,
                        std::string const &columnIndexOrder,
                        std::string const &delimiter = " ");

  template <typename F>
  static void do_write(const std::string &name,
                       const std::string fileName,
                       Tensor<F> *A,
                       const bool binary_p,
                       const std::string rowIndexOrder,
                       const std::string columnIndexOrder,
                       const std::string delimiter);

protected:
  template <typename F = real, typename T = Tensor<F>>
  static T *readBinaryHeader(MPI_File &file, int64_t &offset);
  template <typename F = real, typename T = Tensor<F>>
  static T *readTextHeader(Scanner &scanner);
};
} // namespace sisi4s

#endif
