set(CMAKE_CXX_COMPILER mpiicc)

set(
  OPTIMIZATION_FLAGS
  -xAVX
  -O3
  -ipo
)

set(
  ADDITIONAL_FLAGS
  -DINTEL_COMPILER
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


