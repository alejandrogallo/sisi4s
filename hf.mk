# libint library
LIBINT_PATH = lib/build/${CONFIG}/libint
LIBINT_LIB = -L${LIBINT_PATH}/lib -lint2
LIBINT_INCLUDE = -I${LIBINT_PATH}/include/ -I${LIBINT_PATH}/include/libint2

# eigen library
EIGEN_INCLUDE = -Ilib/build/${CONFIG}/eigen/include/eigen3/

INCLUDE += ${EIGEN_INCLUDE} ${LIBINT_INCLUDE}
LIBS += ${LIBINT_LIB}
