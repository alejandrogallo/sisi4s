include(ExternalProject)

ExternalProject_Add(libint_project
GIT_REPOSITORY    https://github.com/evaleev/libint.git
GIT_TAG           v2.4.2
PREFIX            ${CMAKE_BINARY_DIR}/lib/build/
SOURCE_DIR        ${CMAKE_BINARY_DIR}/lib/src/libint
BINARY_DIR        ${CMAKE_BINARY_DIR}/lib/src/libint
UPDATE_COMMAND    ""
CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/lib/src/libint/autogen.sh &&
                  ${CMAKE_BINARY_DIR}/lib/src/libint/configure --prefix=${CMAKE_BINARY_DIR}/lib/build/libint
BUILD_COMMAND     make
INSTALL_COMMAND   make install
TEST_COMMAND      ""
)

ExternalProject_Get_Property(libint_project SOURCE_DIR)
ExternalProject_Get_Property(libint_project PREFIX)
message("LIBINT SOURCE_DIR = ${SOURCE_DIR}")
set(LIBINT_LIB_DIR ${PREFIX}/libint/lib)
set(LIBINT_INCLUDE_DIR ${PREFIX}/libint/include)
