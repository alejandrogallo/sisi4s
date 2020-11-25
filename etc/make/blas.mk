# BLAS and LAPACK library, expects BLAS_PATH to be set
BLAS_LIB ?= -L${BLAS_PATH}/lib -lopenblas
BLAS_INCLUDE ?=-I${BLAS_PATH}/include

CC4S_INCLUDE     += ${BLAS_INCLUDE}
LIBS_STATIC += ${BLAS_LIB}

