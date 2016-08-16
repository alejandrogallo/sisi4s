# default configuration, override with
#   make all CONFIG=(icc|gcc)
CONFIG=icc

include config.${CONFIG}
include Objects

# goals:
.DEFAULT_GOAL := all
clean:
	rm -rf build/${CONFIG}/bin/*
	rm -rf build/${CONFIG}/obj/*

# primary target
all: build/${CONFIG}/bin/${TARGET}

# generate documentation
doc:
	doxygen

# copy binary to installation directory
install: build/${CONFIG}/bin/${TARGET}
	mkdir -p ${INSTALL}
	cp build/${CONFIG}/bin/${TARGET} ${INSTALL}

# build dependencies only
depend: ${OBJECTS:.o=.d} build/${CONFIG}/obj/${TARGET}.d


# retrieve build environment
VERSION:=$(shell git describe --all --dirty --long)
DATE:=$(shell git log -1 --format="%cd")
COMPILER_VERSION:=$(shell ${CXX} --version | head -n 1)

# add build environment specifics to INCLUDE and to OPTIONS
INCLUDE+=-Isrc
OPTIONS+= -std=c++0x -Wall -fmax-errors=3 -g \
-D_POSIX_C_SOURCE=200112L \
-D__STDC_LIMIT_MACROS -DFTN_UNDERSCORE=1 -DCC4S_VERSION=\"${VERSION}\" \
"-DCC4S_DATE=\"${DATE}\"" \
"-DCOMPILER_VERSION=\"${COMPILER_VERSION}\""


# create a dependency for object file
build/${CONFIG}/obj/%.d: src/%.cxx
	mkdir -p $(dir $@)
	${CXX} -MM ${OPTIONS} ${INCLUDE} -c src/$*.cxx | \
	  sed 's#[^ :]*\.o[ :]*#build/${CONFIG}/obj/$*.o $@: #g' > $@

# keep dependency files
.PRECIOUS: build/${CONFIG}/obj/%.o ${OBJECTS}

# include created dependencies
-include ${OBJECTS:.o=.d}
-include build/${CONFIG}/obj/${TARGET}.d


# compile a object file
build/${CONFIG}/obj/%.o: build/${CONFIG}/obj/%.d
	mkdir -p $(dir $@)
	${CXX} ${OPTIONS} ${OPTIMIZE} ${INCLUDE} -c src/$*.cxx -o $@

# keep object files
.PRECIOUS: build/${CONFIG}/obj/%.o ${OBJECTS}

# compile and link executable
build/${CONFIG}/bin/%: build/${CONFIG}/obj/%.o ${OBJECTS}
	mkdir -p $(dir $@)
	${CXX} ${OPTIONS} ${OPTIMIZE} ${OBJECTS} build/${CONFIG}/obj/${TARGET}.o ${INCLUDE} ${LIBS} -o $@

