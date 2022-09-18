ECL_GIT_REPOSITORY ?= https://gitlab.com/embeddable-common-lisp/ecl
ECL_STATIC_LIB = ${ECL_BUILD_PATH}/libecl.a
ECL_CONFIGURE = $(abspath $(ECL_SRC_PATH))/src/configure
ECL_MAKEFILE = $(abspath $(ECL_BUILD_PATH))/Makefile

$(ECL_CONFIGURE):
	mkdir -p $(ECL_SRC_PATH)
	git clone --depth=1 -b ${ECL_COMMIT} $(ECL_GIT_REPOSITORY) $(ECL_SRC_PATH)

${ECL_MAKEFILE}: $(ECL_CONFIGURE)
	mkdir -p $(dir $@)
	cd ${ECL_BUILD_PATH} && \
	$< \
		--prefix=$(abspath ${ECL_BUILD_PATH}) \
		CXX="$(CXX)"     \
		--with-asdf      \
		--with-cmp       \
		--disable-soname \
		--disable-shared \
		CXXFLAGS="-std=c++11"

${ECL_STATIC_LIB}: $(ECL_MAKEFILE)
	cd $(<D)
	$(MAKE)

ecl: ${ECL_STATIC_LIB}
.PHONY: ecl
