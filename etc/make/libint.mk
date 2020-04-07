# libint library
LIBINT_PATH    ?= lib/build/${CONFIG}/libint
LIBINT_LIB     ?= -L${LIBINT_PATH}/lib -lint2
LIBINT_INCLUDE ?= -I${LIBINT_PATH}/include/ -I${LIBINT_PATH}/include/libint2

INCLUDE        += ${LIBINT_INCLUDE}
LIBS_STATIC    += ${LIBINT_LIB}
