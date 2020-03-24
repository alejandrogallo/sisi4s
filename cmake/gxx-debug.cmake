set(CMAKE_CXX_COMPILER mpicxx)

include(lapack.cmake)
include(scalapack.cmake)
include(libint.cmake)
include(eigen.cmake)

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
  ctf scalapack lapack blas gfortran int2
)

set(
  CC4S_LINK_DIRECTORIES
  ${SCALAPACK_LIB_DIR}
  ${LAPACK_LIB_DIR}
  ${LIBINT_LIB_DIR}
)

set(
  CC4S_INCLUDE_DIRECTORIES
  ${LIBINT_INCLUDE_DIR}
  ${EIGEN_INCLUDE_DIR}
)

set(
  LINKER_FLAGS "-fopenmp"
)
