HF = yes
SCALAPACK_PATH = /opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/gcc-9.1.0/netlib-scalapack-2.0.2-vhgfz5yny2675ijoyxbep35tmk55a3zj
#LIBINT_PATH = $(PWD)/../libint/build/

# destination path for installation
INSTALL=~/bin/cc4s/${CONFIG}

ifeq ($(HF), yes)
include etc/make/eigen.mk
include etc/make/libint.mk
endif

include etc/make/yaml.mk
include etc/make/ctf.mk
