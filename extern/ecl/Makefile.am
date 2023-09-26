
ECL_BUILD_PATH = $(builddir)/src/$(ECL_COMMIT)/build/
ECL_SRC_PATH = $(builddir)/src/$(ECL_COMMIT)/src/
ECL_INSTALL_PATH = $(abs_builddir)

ECL_GIT_REPOSITORY ?= https://gitlab.com/embeddable-common-lisp/ecl
ECL_STATIC_LIB = ${ECL_BUILD_PATH}/libecl.a
ECL_CONFIGURE = $(abspath $(ECL_SRC_PATH))/src/configure
ECL_MAKEFILE = $(abspath $(ECL_BUILD_PATH))/Makefile

################################################################################
include_HEADERS = include/ecl/ecl.h
lib_LIBRARIES = lib/libecl.a

include/ecl/ecl.h: lib/libecl.a
lib/libecl.a: ${ECL_STATIC_LIB}
	$(MAKE) install -C $(ECL_BUILD_PATH)

################################################################################

$(ECL_CONFIGURE):
	mkdir -p $(ECL_SRC_PATH)
	git clone --depth=1 -b ${ECL_COMMIT} $(ECL_GIT_REPOSITORY) $(ECL_SRC_PATH)

${ECL_MAKEFILE}: $(ECL_CONFIGURE)
	mkdir -p $(dir $@)
	cd ${ECL_BUILD_PATH} && \
	$< \
		--prefix=$(ECL_INSTALL_PATH) \
		CXX="$(CXX)"     \
		--with-asdf      \
		--with-cmp       \
		--with-dffi=included       \
		--with-gmp=included       \
		--disable-soname \
		--disable-shared \
		--enable-c99complex       \
		--enable-bohem=included       \
		CXXFLAGS="-std=c++11"

${ECL_STATIC_LIB}: $(ECL_MAKEFILE)
	$(MAKE) -C $(<D)

ecl: ${ECL_STATIC_LIB}
.PHONY: ecl
