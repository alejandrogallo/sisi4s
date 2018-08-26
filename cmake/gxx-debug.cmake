set(CMAKE_CXX_COMPILER mpicxx)

include(lapack.cmake)
include(scalapack.cmake)

set(
  OPTIMIZATION_FLAGS
  -g
  -O0
)

set(
  ADDITIONAL_FLAGS
  -DDEBUG
)

set(
  LIBS
  ctf scalapack lapack blas gfortran
)

set(
  CC4S_LINK_DIRECTORIES
  ${SCALAPACK_LIB_DIR}
  ${LAPACK_LIB_DIR}
)

set(
  LINKER_FLAGS "-fopenmp"
)
