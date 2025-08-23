BT2_DIR = nvbio

MAKEFLAGS="-j 8"
export MAKEFLAGS

CXX = g++

.PHONY: all clean

all:
	cd "../src/nvbio/build" && rm -f CMakeCache.txt && cmake ..
	(cd "../src/nvbio/build" && $(MAKE))
	cd "../src/nvbio/build" && cp nvBWT/nvBWT ../../../inst
	cd "../src/nvbio/build" && cp nvBowtie/nvBowtie ../../../inst
	cd "../src/nvbio/build" && cp nvbio-test/nvbio-test ../../../inst
	

clean:
	(cd "../src/nvbio/build" && rm -f CMakeCache.txt && ($(MAKE) clean -f Makefile))
	cd "../src/nvbio/build" && rm -f *.so
	cd "../src/nvbio/build" && rm -f *.o
	cd "../src/nvbio/build" && rm -f nvBWT
	cd "../src/nvbio/build" && rm -f nvBowtie
	cd "../src/nvbio/build" && rm -f nvbio-test
	cd "../int" && rm -f nvBWT
	cd "../int" && rm -f nvBowtie
	cd "../int" && rm -f nvbio-test
