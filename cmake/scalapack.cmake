include(ExternalProject)

ExternalProject_Add(scalapack_project
URL       http://www.netlib.org/scalapack/scalapack.tgz
PREFIX            ${CMAKE_BINARY_DIR}/lib/src/
SOURCE_DIR        ${CMAKE_BINARY_DIR}/lib/src/scalapack
BINARY_DIR        ${CMAKE_BINARY_DIR}/lib/build/scalapack
CMAKE_ARGS        ""
BUILD_COMMAND     make
INSTALL_COMMAND   ""
UPDATE_COMMAND   ""
TEST_COMMAND      ""
)

ExternalProject_Get_Property(scalapack_project SOURCE_DIR)
ExternalProject_Get_Property(scalapack_project BINARY_DIR)
message("LAPACK SOURCE_DIR = ${SOURCE_DIR}")
message("LAPACK BINARY_DIR = ${BINARY_DIR}")
set(SCALAPACK_LIB_DIR ${BINARY_DIR}/lib)
