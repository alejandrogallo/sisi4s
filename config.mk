CXX=mpiicc
CXXOPTIONS=-openmp -ipo -DINTEL_COMPILER
CXXOPTIMIZE=-march=core-avx2
LIBS=-mkl
#LIBS=-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm
