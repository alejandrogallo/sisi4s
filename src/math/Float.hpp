#ifndef FLOAT_DEFINED
#define FLOAT_DEFINED

#ifndef INTEL_COMPILER
#include <ostream>
#endif


// TODO: use configuration for setting default float type sizes in bits
#define DEFAULT_FLOAT_BIT_SIZE 64
#define MACHINE_FLOAT_BIT_SIZE 64

namespace cc4s {
  template <int FloatSize>
  class FloatTypes;

  template <>
  class FloatTypes<32> {
  public:
    typedef float type;
  };

  template <>
  class FloatTypes<64> {
  public:
    typedef double type;
  };

  // define explicit size float types
  typedef FloatTypes<32>::type Float32;
  typedef FloatTypes<64>::type Float64;

  // define machine supported float as real type
  typedef FloatTypes<
    MACHINE_FLOAT_BIT_SIZE < DEFAULT_FLOAT_BIT_SIZE ?
      MACHINE_FLOAT_BIT_SIZE : DEFAULT_FLOAT_BIT_SIZE
  >::type real;

}
#endif

