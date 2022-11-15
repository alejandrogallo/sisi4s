EIGEN_BUILD_PATH = $(builddir)/src/$(EIGEN_BRANCH)/build
EIGEN_SRC_PATH   = $(builddir)/src/$(EIGEN_BRANCH)/src
EIGEN_CPPFLAGS   ?= -I$(EIGEN_BUILD_PATH)/include/eigen3
CPPFLAGS         += $(EIGEN_CPPFLAGS)

EIGEN_GIT_REPOSITORY ?= https://gitlab.com/libeigen/eigen.git
EIGEN_CMAKE = $(EIGEN_SRC_PATH)/CMakeLists.txt
EIGEN_MAKE = $(EIGEN_BUILD_PATH)/Makefile

EIGEN_MAIN = $(abs_builddir)/include

EIGEN_INSTALL = cd $(EIGEN_BUILD_PATH) && $(MAKE) install

################################################################################

nobase_include_HEADERS = $(EIGEN_MAIN)
${EIGEN_MAIN}: $(EIGEN_MAKE)
	cd ${<D} && $(MAKE) install

################################################################################

$(EIGEN_CMAKE):
	mkdir -p $(@D)
	git clone --depth=1 -b $(EIGEN_BRANCH) $(EIGEN_GIT_REPOSITORY) $(@D)
	touch $@

${EIGEN_MAKE}: $(EIGEN_CMAKE)
	@echo ${EIGEN_MAIN} aus $(EIGEN_CMAKE)
	mkdir -p $(@D)
	cd $(@D) ; \
	cmake -DCMAKE_INSTALL_PREFIX=$(abs_builddir) $(abspath $(<D))

################################################################################

eigen-install: $(EIGEN_MAKE)

fetch: $(EIGEN_CMAKE)
eigen: $(EIGEN_MAIN)
clean-local:
	-rm -r $(EIGEN_BUILD_PATH)
	-rm -r include share

.PHONY: all eigen-clean eigen fetch eigen-install
