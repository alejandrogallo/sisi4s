LIBINT_MAX_AM ?= 4
LIBINT_GIT_REPOSITORY ?= https://github.com/evaleev/libint
LIBINT_STATIC_LIB = ${LIBINT_BUILD_PATH}/lib/libint2.la
LIBINT_AUTOGEN = $(LIBINT_SRC_PATH)/autogen.sh
LIBINT_CONFIGURE = $(abspath $(LIBINT_SRC_PATH))/configure

$(LIBINT_AUTOGEN):
	git clone --depth=1 -b ${LIBINT_COMMIT} $(LIBINT_GIT_REPOSITORY) $(@D)

$(LIBINT_CONFIGURE): $(LIBINT_AUTOGEN)
	cd $(@D) && ./autogen.sh

${LIBINT_STATIC_LIB}: $(LIBINT_CONFIGURE)
	mkdir -p $(dir $@)
	cd ${LIBINT_BUILD_PATH} && \
	$< \
		--with-max-am=${LIBINT_MAX_AM} \
		--prefix=$(abspath ${LIBINT_BUILD_PATH}) \
		--with-boost=${BOOST_ROOT} \
		CXX="$(CXX)" \
		CXXFLAGS="-std=c++11" && \
	$(MAKE) install

libint: ${LIBINT_STATIC_LIB}
.PHONY: libint
