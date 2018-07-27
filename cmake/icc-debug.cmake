set(CMAKE_CXX_COMPILER mpiicc)

set(
  OPTIMIZATION_FLAGS
  -DDEBUG
  -O0
  -g
)

set(
  ADDITIONAL_FLAGS
  -DINTEL_COMPILER
  -Qoption,cpp,--extended_float_types
  -mkl
  -openmp
)

set(
  LIBS
  ctf mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64
)

set(
  LINKER_FLAGS "-mkl -openmp"
)


