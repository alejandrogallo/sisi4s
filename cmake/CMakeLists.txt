cmake_minimum_required(VERSION 3.5) # setting this is required
project(Cc4s)            # this sets the project name


set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
set(CONFIG gxx-debug)

# Check for openmp and add the openmp flag also to the linker
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Find lapack libraries
find_path (
  SCALAPACK_DIR libscalapack.a
  HINTS ENV SCALAPACK_DIR
  PATHS "/usr/lib/scalapack" "/usr/local/lib/scalapack"
  DOC "ScaLapack Directory"
)


set(CTF "${PROJECT_SOURCE_DIR}/lib/build/${CONFIG}/ctf")
find_path (
  CTF_DIR lib/libctf.a
  HINTS ENV CTF_DIR
  PATHS ${CTF} "/usr/lib/ctf" "/usr/local/lib/ctf"
  DOC "CTF Directory"
)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${CTF_DIR}/include
)

link_directories(
  ${SCALAPACK_DIR} ${CTF_DIR}/lib /usr/lib
)

set(CMAKE_CXX_COMPILER mpicxx)
string(TIMESTAMP CC4S_DATE)
execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} --version
  COMMAND head -n1
  COMMAND tr -d "\n"
  OUTPUT_VARIABLE COMPILER_VERSION
)
execute_process(
  COMMAND git describe --all --dirty --long
  COMMAND tr -d "\n"
  OUTPUT_VARIABLE VERSION
)


message("CONFIG           = ${CONFIG}")
message("CTF_DIR          = ${CTF_DIR}")
message("SCALAPACK_DIR    = ${SCALAPACK_DIR}")
message("CC4S_DATE        = ${CC4S_DATE}")
message("VERSION          = ${VERSION}")
message("CXX              = ${CMAKE_CXX_COMPILER}")
message("COMPILER_VERSION = ${COMPILER_VERSION}")
message("OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")

# DEFINE THE SOURCES FOR THE CC4S BINARY
file(
  GLOB
  Cc4s_sources
  ${PROJECT_SOURCE_DIR}/src/*
  ${PROJECT_SOURCE_DIR}/src/algorithms/*
  ${PROJECT_SOURCE_DIR}/src/extern/*
  ${PROJECT_SOURCE_DIR}/src/gte/*
  ${PROJECT_SOURCE_DIR}/src/math/*
  ${PROJECT_SOURCE_DIR}/src/mixers/*
  ${PROJECT_SOURCE_DIR}/src/tcc/*
  ${PROJECT_SOURCE_DIR}/src/util/*
)


# DECLARE AN EXECUTABLE CALLED Cc4s
add_executable(Cc4s ${Cc4s_sources})

# declare libraries needed for executable
target_link_libraries(
  Cc4s
  ctf scalapack reflapack refblas gfortran
)

# set compile flags for executable Cc4s
target_compile_options(
  Cc4s
  PUBLIC
  ${OpenMP_CXX_FLAGS}
  -fmax-errors=3
  -Wall
  -std=c++11
  -g
  -O0
  -DDEBUG
  -D_POSIX_C_SOURCE=200112L
  -D__STDC_LIMIT_MACROS
  -DFTN_UNDERSCORE=1
  -DCC4S_VERSION="${VERSION}"
  -DCC4S_DATE="${CC4S_DATE}"
  -DCOMPILER_VERSION="${COMPILER_VERSION}"
)

set_target_properties(
  Cc4s
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
