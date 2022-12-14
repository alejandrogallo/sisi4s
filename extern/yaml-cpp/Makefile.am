AUTOMAKE_OPTIONS = subdir-objects

YAML_BUILD_PATH     = $(abs_builddir)/src/$(YAML_COMMIT)/build
YAML_SRC_PATH       = $(abs_builddir)/src/$(YAML_COMMIT)/src

YAML_GIT_REPOSITORY ?= https://github.com/cc4s/yaml-cpp.git
YAML_STATIC_LIB = $(YAML_BUILD_PATH)/libyaml-cpp.a
YAML_CMAKE = $(YAML_SRC_PATH)/CMakeLists.txt
YAML_MAKE = $(YAML_BUILD_PATH)/Makefile
CMAKE ?= cmake

YAML_INSTALL = cd $(YAML_BUILD_PATH) && $(MAKE) install

################################################################################

lib_LIBRARIES = lib/libyaml-cpp.a
lib_libyaml_cpp_a_SOURCES = $(YAML_STATIC_LIB)
lib/libyaml-cpp.a: $(YAML_STATIC_LIB)
	mkdir -p lib
	$(YAML_INSTALL)
	cp $< $@

################################################################################

$(YAML_CMAKE):
	mkdir -p $(@D)
	git clone  ${YAML_GIT_REPOSITORY} $(@D)
	cd $(@D) && git checkout $(YAML_COMMIT)

$(YAML_MAKE): $(YAML_CMAKE)
	mkdir -p $(YAML_BUILD_PATH)
	cd $(YAML_BUILD_PATH) && \
	$(CMAKE) \
		-DCMAKE_INSTALL_PREFIX=$(abs_builddir) \
		-DYAML_CPP_BUILD_TESTS:BOOL=OFF \
		-DYAML_CPP_BUILD_TOOLS:BOOL=OFF \
		-DCMAKE_CXX_FLAGS:STRING="-fPIE" \
		-DCMAKE_CXX_COMPILER:STRING="$(MPICXX)" \
		$(abspath $(YAML_SRC_PATH))

$(YAML_STATIC_LIB): $(YAML_MAKE)
	$(info Compiling $@)
	cd $(@D) && $(MAKE)

################################################################################

yaml-cpp: $(YAML_STATIC_LIB)

clean-local:
	rm -rf $(YAML_BUILD_PATH)
	-rm -rf include
	-rm -rf lib
	-rm -rf share
	-rm -rf lib64

fetch: $(YAML_CMAKE)

.PHONY: yaml-cpp fetch
