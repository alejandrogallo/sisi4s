set(CMAKE_CXX_COMPILER mpicxx)

set(
  OPTIMIZATION_FLAGS
  -march=native
  -flto
  -Ofast
)

set(
  ADDITIONAL_FLAGS
)

set(
  LIBS
  ctf scalapack reflapack refblas gfortran
)

# Find lapack libraries
find_path (
  SCALAPACK_DIR libscalapack.a
  HINTS ENV SCALAPACK_DIR
  PATHS "/usr/lib/scalapack" "/usr/local/lib/scalapack"
  DOC "ScaLapack Directory"
)

set(
  LINKER_FLAGS "-fopenmp"
)
