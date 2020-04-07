CONFIG ?= gxx
include config.${CONFIG}

HF ?= yes
THERMAL ?= no

# destination path for installation
INSTALL ?= ~/bin/cc4s/${CONFIG}
TARGET ?= Cc4s

ifeq ($(HF), yes)
include etc/make/eigen.mk
include etc/make/libint.mk
endif

include etc/make/yaml.mk
include etc/make/ctf.mk

ifdef INTEL_COMPILER
LIBS += ${LIBS_STATIC} ${LIBS_DYNAMIC}
else
LIBS += -Wl,-Bstatic ${LIBS_STATIC} -Wl,-Bdynamic ${LIBS_DYNAMIC}
endif
