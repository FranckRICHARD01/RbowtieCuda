BT2_DIR = nvbio

MAKEFLAGS="-j 8 -l 8"
export MAKEFLAGS

CXX = g++

BUILD_DIR = ../src/nvbio/build

.PHONY: all clean

all:
	@mkdir -p "$(BUILD_DIR)"
	cd "$(BUILD_DIR)" && rm -f CMakeCache.txt && cmake ..
	(cd "$(BUILD_DIR)" && $(MAKE))
	cd "$(BUILD_DIR)" && cp nvBWT/nvBWT ../../../inst
	cd "$(BUILD_DIR)" && cp nvBowtie/nvBowtie ../../../inst
	cd "$(BUILD_DIR)" && cp nvbio-test/nvbio-test ../../../inst

clean:
	(cd "$(BUILD_DIR)" && rm -f CMakeCache.txt && ($(MAKE) clean -f Makefile))
	cd "$(BUILD_DIR)" && rm -f *.so
	cd "$(BUILD_DIR)" && rm -f *.o
	cd "$(BUILD_DIR)" && rm -f nvBWT nvBowtie nvbio-test
	cd ../../../inst && rm -f nvBWT nvBowtie nvbio-test
