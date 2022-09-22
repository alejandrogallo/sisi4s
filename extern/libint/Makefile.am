LIBINT_MAX_AM ?= 4
LIBINT_GIT_REPOSITORY ?= https://github.com/evaleev/libint
LIBINT_STATIC_LIB = ${LIBINT_BUILD_PATH}/lib/.libs/libint2.a
LIBINT_AUTOGEN = $(LIBINT_SRC_PATH)/autogen.sh
LIBINT_CONFIGURE = $(abspath $(LIBINT_SRC_PATH))/configure

LIBINT_BUILD_PATH = $(builddir)/src/$(LIBINT_COMMIT)/build
LIBINT_SRC_PATH   = $(builddir)/src/$(LIBINT_COMMIT)/src

LIBINT_INSTALL = cd ${LIBINT_BUILD_PATH} && $(MAKE) install

################################################################################

lib_LIBRARIES = lib/libint2.a
lib/libint2.a: $(LIBINT_STATIC_LIB)
	$(LIBINT_INSTALL)

nobase_include_HEADERS = $(patsubst $(LIBINT_BUILD_PATH)/%,\
                             ./%,\
                             $(shell find $(LIBINT_BUILD_PATH)/include -type f))
$(nobase_include_HEADERS):
	$(LIBINT_INSTALL)


################################################################################

$(LIBINT_AUTOGEN):
	git clone --depth=1 -b ${LIBINT_COMMIT} $(LIBINT_GIT_REPOSITORY) $(@D)

$(LIBINT_CONFIGURE): $(LIBINT_AUTOGEN)
	cd $(@D) && ./autogen.sh

${LIBINT_STATIC_LIB}: $(LIBINT_CONFIGURE)
	mkdir -p $(dir $@)
	cd ${LIBINT_BUILD_PATH} && \
	$< \
		--with-max-am=${LIBINT_MAX_AM} \
		--prefix=$(abs_builddir) \
		--with-boost=${BOOST_ROOT} \
		CXX="$(CXX)" \
		CXXFLAGS="-std=c++11" && \
	$(MAKE)

################################################################################


clean-local:
	-rm -rf include
	-rm -rf lib
	-rm -rf share

libint: ${LIBINT_STATIC_LIB}
fetch: $(LIBINT_AUTOGEN)

.PHONY: libint fetch
