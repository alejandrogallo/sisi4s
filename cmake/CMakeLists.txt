cmake_minimum_required(VERSION 2.8) # setting this is required
project(Cc4s CXX)            # this sets the project name

set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
if("${CONFIG}" STREQUAL "")
  set(CONFIG gxx-debug)
endif()

include(${CONFIG}.cmake)
include(ctf.cmake)
include(libint.cmake)

# Check for openmp and add the openmp flag also to the linker
find_package(OpenMP REQUIRED)

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

message("CONFIG             = ${CONFIG}")
message("CTF_DIR            = ${CTF_DIR}")
message("CC4S_DATE          = ${CC4S_DATE}")
message("VERSION            = ${VERSION}")
message("CXX                = ${CMAKE_CXX_COMPILER}")
message("COMPILER_VERSION   = ${COMPILER_VERSION}")
message("OpenMP_CXX_FLAGS   = ${OpenMP_CXX_FLAGS}")
message("OPTIMIZATION_FLAGS = ${OPTIMIZATION_FLAGS}")
message("ADDITIONAL_FLAGS   = ${ADDITIONAL_FLAGS}")

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

# directories to look for libraries
link_directories(
  ${CC4S_LINK_DIRECTORIES}
  ${CTF_DIR}/lib /usr/lib
)

# declare an executable called Cc4s
add_executable(Cc4s ${Cc4s_sources})

target_include_directories(
  Cc4s
  PUBLIC
  ${PROJECT_SOURCE_DIR}/src
  ${CTF_DIR}/include
)

# declare libraries needed for executable
target_link_libraries(Cc4s ${LIBS})

# set compile flags for executable Cc4s
target_compile_options(
  Cc4s
  PUBLIC
  ${OpenMP_CXX_FLAGS}
  -fmax-errors=3
  -Wall
  -std=c++11
  ${OPTIMIZATION_FLAGS}
  ${ADDITIONAL_FLAGS}
  -D_POSIX_C_SOURCE=200112L
  -D__STDC_LIMIT_MACROS
  -DFTN_UNDERSCORE=1
  -DCC4S_VERSION="${VERSION}"
  -DCC4S_DATE="${CC4S_DATE}"
  -DCOMPILER_VERSION="${COMPILER_VERSION}"
)

set_target_properties(
  Cc4s
  PROPERTIES LINK_FLAGS "${LINKER_FLAGS}"
)
set_target_properties(
  Cc4s
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
