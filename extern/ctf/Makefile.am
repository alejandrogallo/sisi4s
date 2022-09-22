CTF_GIT_REPOSITORY ?= https://gitlab.cc4s.org/cc4s/ctf.git

CTF_BUILD_PATH     ?= $(builddir)/src/$(CTF_COMMIT)/build
CTF_SRC_PATH       ?= $(builddir)/src/$(CTF_COMMIT)/src

CTF_CONFIG_FLAGS = CXX="$(MPICXX) -std=c++11" \
                   LIB_PATH="$(LDFLAGS)" \
                   LIBS="${LIBS}" \
                   CXXFLAGS="$(CXXFLAGS)" \
                   --no-dynamic

CTF_STATIC_LIB = $(CTF_BUILD_PATH)/lib/libctf.a
CTF_SHARED_LIB = $(CTF_BUILD_PATH)/lib/libctf.so
CTF_GIT_REPOSITORY ?= https://github.com/cyclops-community/ctf
CTF_CONFIGURE = $(CTF_SRC_PATH)/configure
## CTF_INSTALL = cd $(CTF_BUILD_PATH) && $(MAKE) install

################################################################################

lib_LIBRARIES = lib/libctf.a
lib/libctf.a: $(CTF_STATIC_LIB)
	mkdir -p lib
	cp $(CTF_BUILD_PATH)/lib/libctf.a lib/

include_HEADERS = include/ctf.hpp
include/ctf.hpp: $(CTF_STATIC_LIB)
	mkdir -p include
	cp $(CTF_BUILD_PATH)/include/ctf.hpp $@

################################################################################

$(CTF_CONFIGURE):
	mkdir -p $(@D)
	git clone $(CTF_GIT_REPOSITORY) $(@D)
	cd $(@D) && git checkout $(CTF_COMMIT)

$(CTF_BUILD_PATH)/Makefile: $(CTF_SRC_PATH)/configure
	mkdir -p $(CTF_BUILD_PATH)
	cd $(CTF_BUILD_PATH) && \
		$(abspath $(CTF_SRC_PATH))/configure \
											--install-dir=$(abs_builddir) \
											$(CTF_CONFIG_FLAGS)

$(CTF_STATIC_LIB): $(CTF_BUILD_PATH)/Makefile
	$(info Compiling $@)
	cd $(CTF_BUILD_PATH) && $(MAKE) ctflib

################################################################################

clean-local:
	-rm -rf include
	-rm -rf lib

ctf: $(CTF_STATIC_LIB)

fetch: $(CTF_CONFIGURE)

.PHONY: ctf fetch
