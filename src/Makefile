BT2_DIR = nvbio

MAKEFLAGS="-j 4"

CXX = g++

.PHONY: all clean

all:
	mkdir -p ../inst
	(cd "../build" && rm -f CMakeCache.txt && cmake "../src/"$(BT2_DIR) .)
	(cd "../build" && $(MAKE))
	cd "../build" && cp nvBWT/nvBWT ../inst
	cd "../build" && cp nvBowtie/nvBowtie ../inst
	cd "../build" && cp nvbio-test/nvbio-test ../inst
	

clean:
	(cd "../build" && rm -f CMakeCache.txt && ($(MAKE) clean -f Makefile))
	rm -f *.so
	rm -f *.o
	rm -f nvBWT
	rm -f nvBowtie
	rm -f nvbio-test
	rm -f ../inst/nvBWT
	rm -f ../inst/nvBowtie
	rm -f ../inst/nvbio-test
