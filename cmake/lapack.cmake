include(ExternalProject)

ExternalProject_Add(lapack_project
URL       http://www.netlib.org/lapack/lapack-3.8.0.tar.gz
PREFIX            ${CMAKE_BINARY_DIR}/lib/src/
SOURCE_DIR        ${CMAKE_BINARY_DIR}/lib/src/lapack
BINARY_DIR        ${CMAKE_BINARY_DIR}/lib/build/lapack
CMAKE_ARGS        ""
BUILD_COMMAND     make
INSTALL_COMMAND   ""
UPDATE_COMMAND   ""
TEST_COMMAND      ""
)

ExternalProject_Get_Property(lapack_project SOURCE_DIR)
ExternalProject_Get_Property(lapack_project BINARY_DIR)
message("LAPACK SOURCE_DIR = ${SOURCE_DIR}")
message("LAPACK BINARY_DIR = ${BINARY_DIR}")
set(LAPACK_LIB_DIR ${BINARY_DIR}/lib)

