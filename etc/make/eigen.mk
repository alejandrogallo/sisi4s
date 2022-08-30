EIGEN_GIT_REPOSITORY ?= https://gitlab.com/libeigen/eigen.git
EIGEN_CMAKE = $(EIGEN_SRC_PATH)/CMakeLists.txt
EIGEN_MAKE = $(EIGEN_BUILD_PATH)/Makefile

EIGEN_MAIN = $(EIGEN_BUILD_PATH)/include/eigen3/Eigen/Eigen

$(EIGEN_CMAKE):
	mkdir -p $(@D)
	git clone --depth=1 -b $(EIGEN_BRANCH) $(EIGEN_GIT_REPOSITORY) $(@D)
	touch $@

${EIGEN_MAKE}: $(EIGEN_CMAKE)
	@echo ${EIGEN_MAIN} aus $(EIGEN_CMAKE)
	mkdir -p $(@D)
	cd $(@D) ; \
	cmake -DCMAKE_INSTALL_PREFIX=$(abspath $(@D)) $(abspath $(<D))

${EIGEN_MAIN}: $(EIGEN_MAKE)
	cd ${<D} && $(MAKE) install
	test -f $(@) && touch $@

eigen: $(EIGEN_MAIN)
eigen-clean:
	rm -r $(EIGEN_BUILD_PATH)

.PHONY: eigen
