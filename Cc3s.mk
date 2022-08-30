EXTERNAL_DEPENDENCIES =
extern: $(EXTERNAL_DEPENDENCIES)
.PHONY: extern

if WITH_BUILD_CTF
# Cyclops Tensor Framework ====================================================
CTF_COMMIT         ?= 53ae5daad851bf3b198ebe1fa761c13b12291116
CTF_BUILD_PATH     ?= $(top_builddir)/extern/build/ctf/$(CTF_COMMIT)
CTF_SRC_PATH       ?= $(top_builddir)/extern/src/ctf/$(CTF_COMMIT)
CTF_CPPFLAGS       ?= -I$(CTF_BUILD_PATH)/include
CTF_LDFLAGS        ?= -L${CTF_BUILD_PATH}/lib -lctf
CTF_GIT_REPOSITORY ?= https://gitlab.cc4s.org/cc4s/ctf.git
EXTERNAL_DEPENDENCIES += ctf
include $(top_srcdir)/etc/make/ctf.mk
endif

if WITH_BUILD_YAML
# yaml-cpp ====================================================================
YAML_COMMIT         ?= c9460110e072df84b7dee3eb651f2ec5df75fb18
YAML_BUILD_PATH     ?= $(top_builddir)/extern/build/yaml-cpp/$(YAML_COMMIT)
YAML_SRC_PATH       ?= $(top_builddir)/extern/src/yaml-cpp/$(YAML_COMMIT)
YAML_LDFLAGS        ?= -L${YAML_BUILD_PATH} -lyaml-cpp
YAML_CPPFLAGS       ?= -I${YAML_BUILD_PATH}/include
YAML_GIT_REPOSITORY ?= https://gitlab.cc4s.org/cc4s/yaml-cpp.git
EXTERNAL_DEPENDENCIES += yamlcpp
include $(top_srcdir)/etc/make/yaml.mk
endif

if WITH_BUILD_LIBINT
LIBINT_COMMIT     ?= v2.6.0
LIBINT_MAX_AM     ?= 4
LIBINT_BUILD_PATH ?= $(top_builddir)/extern/build/libint/$(LIBINT_COMMIT)
LIBINT_SRC_PATH   ?= $(top_builddir)/extern/src/libint/$(LIBINT_COMMIT)
LIBINT_LIB        ?= ${LIBINT_BUILD_PATH}/lib/libint2.a
LIBINT_CPPFLAGS   ?= -I${LIBINT_BUILD_PATH}/include
EXTERNAL_DEPENDENCIES += libint
include $(top_srcdir)/etc/make/libint.mk
endif

if WITH_BUILD_EIGEN
EIGEN_BUILD_PATH ?= $(top_builddir)/extern/build/eigen/$(EIGEN_BRANCH)
EIGEN_SRC_PATH ?= $(top_builddir)/extern/src/eigen/$(EIGEN_BRANCH)
EIGEN_LIB ?= $(EIGEN_BUILD_PATH)/include/eigen3/Eigen/Eigen
EIGEN_CPPFLAGS ?= -I$(EIGEN_BUILD_PATH)/include/eigen3
EIGEN_BRANCH ?= 3.3.7
EXTERNAL_DEPENDENCIES += eigen
include $(top_srcdir)/etc/make/eigen.mk
endif

# default CXXFLAGS
CXXFLAGS  +=                  \
-D_POSIX_C_SOURCE=200112L     \
-D__STDC_LIMIT_MACROS         \
-DFTN_UNDERSCORE=1            \
-DCC4S_VERSION=\"cc3s\" \
"-DCC4S_DATE=\"${DATE}\""     \
"-DCOMPILER_VERSION=\"${COMPILER_VERSION}\"" \
"-DCC4S_MPI_VERSION=\"${CC4S_MPI_VERSION}\""
