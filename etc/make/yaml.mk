# yaml-cpp library
YAML_PATH    ?= lib/build/${CONFIG}/yaml-cpp
YAML_LIB     ?= -L${YAML_PATH} -lyaml-cpp
YAML_INCLUDE ?= -I${YAML_PATH}/include

INCLUDE += ${YAML_INCLUDE}
LIBS_STATIC    += ${YAML_LIB}
