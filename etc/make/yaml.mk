YAML_STATIC_LIB = $(YAML_BUILD_PATH)/libyaml-cpp.a
YAML_CMAKE = $(YAML_SRC_PATH)/CMakeLists.txt
YAML_GIT_REPOSITORY ?= https://github.com/jbeder/yaml-cpp.git
CMAKE ?= cmake

$(YAML_CMAKE):
	mkdir -p $(@D)
	git clone  ${YAML_GIT_REPOSITORY} $(@D)
	cd $(@D) && git checkout $(YAML_COMMIT)

$(YAML_BUILD_PATH)/Makefile: $(YAML_CMAKE)
	mkdir -p $(YAML_BUILD_PATH)
	cd $(YAML_BUILD_PATH) && \
	$(CMAKE) \
		-DCMAKE_INSTALL_PREFIX=$(abspath $(YAML_BUILD_PATH)) \
		-DYAML_CPP_BUILD_TESTS:BOOL=OFF \
		-DYAML_CPP_BUILD_TOOLS:BOOL=OFF \
		$(abspath $(YAML_SRC_PATH))


$(YAML_STATIC_LIB): $(YAML_BUILD_PATH)/Makefile
	$(info Compiling $@)
	cd $(@D) && $(MAKE) install

.PHONY: yamlcpp yamlcpp-clean
yamlcpp: $(YAML_STATIC_LIB)

yamlcpp-clean:
	rm -rf $(YAML_BUILD_PATH)
