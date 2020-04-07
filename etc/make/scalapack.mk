# ScaLAPACK libarary, expects SCALAPACK_PATH to be set
SCALAPACK_LIB ?= -L${SCALAPACK_PATH}/lib -lscalapack

LIBS_DYNAMIC += ${SCALAPACK_LIB}

