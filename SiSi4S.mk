extern: $(EXTERNAL_DEPENDENCIES)
.PHONY: extern

if WITH_BUILD_CTF
# Cyclops Tensor Framework ====================================================
CTF_BUILD_PATH     ?= $(EXTERN_BUILD_PATH)/ctf/$(CTF_COMMIT)
CTF_SRC_PATH       ?= $(EXTERN_SRC_PATH)/ctf/$(CTF_COMMIT)
CTF_CPATH          ?= $(CTF_BUILD_PATH)/include
CTF_CPPFLAGS       ?= -I$(CTF_CPATH)
CTF_LDFLAGS        ?= -L${CTF_BUILD_PATH}/lib -lctf
CTF_GIT_REPOSITORY ?= https://gitlab.cc4s.org/cc4s/ctf.git
CPPFLAGS           += $(CTF_CPPFLAGS)
include $(top_srcdir)/etc/make/ctf.mk
endif

if WITH_BUILD_YAML
# yaml-cpp ====================================================================
YAML_BUILD_PATH     ?= $(EXTERN_BUILD_PATH)/yaml-cpp/$(YAML_COMMIT)
YAML_SRC_PATH       ?= $(EXTERN_SRC_PATH)/yaml-cpp/$(YAML_COMMIT)
YAML_CPATH          ?= ${YAML_BUILD_PATH}/include
YAML_CPPFLAGS       ?= -I${YAML_CPATH}
YAML_LDFLAGS        ?= ${YAML_BUILD_PATH}/libyaml-cpp.a
YAML_GIT_REPOSITORY ?= https://gitlab.cc4s.org/cc4s/yaml-cpp.git
CPPFLAGS            += $(YAML_CPPFLAGS)
include $(top_srcdir)/etc/make/yaml.mk
endif

if WITH_BUILD_LIBINT
LIBINT_MAX_AM     ?= 4
LIBINT_BUILD_PATH ?= $(EXTERN_BUILD_PATH)/libint/$(LIBINT_COMMIT)
LIBINT_SRC_PATH   ?= $(EXTERN_SRC_PATH)/libint/$(LIBINT_COMMIT)
LIBINT_LIB        ?= ${LIBINT_BUILD_PATH}/lib/libint2.a
LIBINT_CPATH      ?= ${LIBINT_BUILD_PATH}/include
LIBINT_CPPFLAGS   ?= -I${LIBINT_CPATH}
CPPFLAGS          += $(LIBINT_CPPFLAGS)
include $(top_srcdir)/etc/make/libint.mk
endif

if WITH_BUILD_EIGEN
EIGEN_BUILD_PATH ?= $(EXTERN_BUILD_PATH)/eigen/$(EIGEN_BRANCH)
EIGEN_SRC_PATH   ?= $(EXTERN_SRC_PATH)/eigen/$(EIGEN_BRANCH)
EIGEN_LIB        ?= $(EIGEN_BUILD_PATH)/include/eigen3/Eigen/Eigen
EIGEN_CPPFLAGS   ?= -I$(EIGEN_BUILD_PATH)/include/eigen3
CPPFLAGS         += $(EIGEN_CPPFLAGS)
include $(top_srcdir)/etc/make/eigen.mk
endif

if WITH_BUILD_ECL
ECL_BUILD_PATH ?= $(EXTERN_BUILD_PATH)/ecl/$(ECL_COMMIT)
ECL_SRC_PATH   ?= $(EXTERN_SRC_PATH)/ecl/$(ECL_COMMIT)
ECL_LIB        ?= ${ECL_BUILD_PATH}/libecl.a
ECL_CPATH      ?= ${ECL_BUILD_PATH}/
ECL_CPPFLAGS   ?= -I${ECL_CPATH}
CPPFLAGS       += $(ECL_CPPFLAGS)
include $(top_srcdir)/etc/make/ecl.mk
endif
