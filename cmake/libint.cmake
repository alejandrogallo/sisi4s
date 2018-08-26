include(ExternalProject)

ExternalProject_Add(libint_project
GIT_REPOSITORY    https://github.com/evaleev/libint
GIT_TAG           master
PREFIX            ${CMAKE_BINARY_DIR}/lib/src/
SOURCE_DIR        ${CMAKE_BINARY_DIR}/lib/src/libint
BINARY_DIR        ${CMAKE_BINARY_DIR}/lib/src/libint
UPDATE_COMMAND    ""
CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/lib/src/libint/autogen.sh &&
                  ${CMAKE_BINARY_DIR}/lib/src/libint/configure
BUILD_COMMAND     make
INSTALL_COMMAND   ""
TEST_COMMAND      ""
)

ExternalProject_Get_Property(libint_project SOURCE_DIR)
ExternalProject_Get_Property(libint_project BINARY_DIR)
message("LIBINT SOURCE_DIR = ${SOURCE_DIR}")
message("LIBINT BINARY_DIR = ${BINARY_DIR}")
set(LIBINT_LIB_DIR ${BINARY_DIR}/lib)
