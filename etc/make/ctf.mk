# Cyclops Tensor Framework libraray
CTF_PATH    ?= lib/build/${CONFIG}/ctf
CTF_LIB     ?= -L${CTF_PATH}/lib -lctf
CTF_INCLUDE ?= -I${CTF_PATH}/include

CC4S_INCLUDE += ${CTF_INCLUDE}
#LIBS    += -Wl,-Bstatic ${CTF_LIB}
LIBS_STATIC  += ${CTF_LIB}
