CONFIG ?= gxx
include config.${CONFIG}

HF ?= yes
THERMAL ?= no

# destination path for installation
INSTALL = ~/bin/cc4s/${CONFIG}

ifeq ($(HF), yes)
include etc/make/eigen.mk
include etc/make/libint.mk
endif

include etc/make/yaml.mk
include etc/make/ctf.mk
include etc/make/blas.mk
include etc/make/scalapack.mk

LIBS += -Wl,-Bstatic ${LIBS_STATIC} -Wl,-Bdynamic ${LIBS_DYNAMIC}
