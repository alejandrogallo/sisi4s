# Cyclops Tensor Framework libraray
CTF_PATH    ?= lib/build/${CONFIG}/ctf
CTF_LIB     ?= -L${CTF_PATH}/lib -lctf
CTF_INCLUDE ?= -I${CTF_PATH}/include

INCLUDE += ${CTF_INCLUDE}
#LIBS    += -Wl,-Bstatic ${CTF_LIB}
LIBS    += ${CTF_LIB}

