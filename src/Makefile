BT2_DIR = nvbio

CXX = g++

BUILD_DIR = ../src/nvbio/build

# ---------------------------------------------------------------------------
# Build parallelism
#
# Two things have to be respected at once:
#
#  * the caller's choice: R 4.6 exports MAKEFLAGS with a job count and a
#    jobserver (" -j16 --jobserver-auth=fifo:/tmp/GMfifo32"), and
#    Bioconductor/CRAN want packages to honour it. Assigning MAKEFLAGS here is
#    forbidden (that is what the old "MAKEFLAGS=-j 8 -l 8" did: it overrode the
#    builder in both directions).
#  * the memory budget: one nvcc job on the thrust/cub-heavy nvbio translation
#    units costs a few GB of RSS, which is what OOM-kills a 16 GB / 4 vCPU
#    GitHub-hosted runner (r-universe). R's own -j count is derived from the
#    CPU count, not from memory, so inheriting it blindly is exactly how the
#    runner gets killed.
#
# So: never assign MAKEFLAGS, and pass an explicit -j to the nvbio sub-make
# equal to min(caller's job count, nproc, available memory / NVBIO_MEM_PER_JOB_KB).
# An explicit -j on a sub-make line bounds the concurrency even when a jobserver
# was inherited (verified with GNU make 4.4.1), so this only ever lowers the
# caller's parallelism, never raises it.
#
# NVBIO_MEM_PER_JOB_KB is a conservative estimate, not a measurement: adjust it
# after timing a real build (see INSTALL). Override the whole thing by hand with
# NVBIO_JOBS=16 R CMD INSTALL ...
# ---------------------------------------------------------------------------
NVBIO_MEM_PER_JOB_KB ?= 3000000
NVBIO_MEMINFO ?= /proc/meminfo

# Job count asked for by the caller, if any ("--jobs=8" form is normalised too).
NVBIO_CALLER_J := $(patsubst --jobs=%,-j%,$(filter -j% --jobs=%,$(MAKEFLAGS)))
# A bare "-j" (unlimited) carries no number to compare against.
NVBIO_CALLER_N := $(patsubst -j%,%,$(filter -j%,$(NVBIO_CALLER_J)))
ifeq ($(strip $(NVBIO_CALLER_N)),)
NVBIO_CALLER_N :=
endif

# Ceiling from RAM and CPU count. Unreadable /proc/meminfo falls back to one
# job, which is slow but never fatal.
NVBIO_MEM_JOBS := $(shell \
  mem=$$(awk '/^MemAvailable/ {print $$2}' $(NVBIO_MEMINFO) 2>/dev/null); \
  cpu=$$(nproc 2>/dev/null || echo 1); \
  j=$$(( $${mem:-0} / $(NVBIO_MEM_PER_JOB_KB) )); \
  [ "$$j" -lt 1 ] && j=1; \
  [ "$$j" -gt "$$cpu" ] && j=$$cpu; \
  echo $$j)

NVBIO_JOBS ?= $(shell \
  j=$(NVBIO_MEM_JOBS); \
  if [ -n "$(NVBIO_CALLER_N)" ] && [ "$(NVBIO_CALLER_N)" -eq "$(NVBIO_CALLER_N)" ] 2>/dev/null \
     && [ "$$j" -gt "$(NVBIO_CALLER_N)" ]; then j=$(NVBIO_CALLER_N); fi; \
  [ "$$j" -lt 1 ] && j=1; \
  echo $$j)

.PHONY: all clean

all:
	@echo "nvbio build: $(NVBIO_JOBS) parallel job(s) ceiling (memory: $(NVBIO_MEM_JOBS), caller: $(or $(NVBIO_CALLER_J),none), per job: $(NVBIO_MEM_PER_JOB_KB) kB)"
	@mkdir -p "$(BUILD_DIR)"
	cd "$(BUILD_DIR)" && rm -f CMakeCache.txt && cmake ..
	(cd "$(BUILD_DIR)" && $(MAKE) -j$(NVBIO_JOBS))
	cd "$(BUILD_DIR)" && cp nvBWT/nvBWT ../../../inst
	cd "$(BUILD_DIR)" && cp nvBowtie/nvBowtie ../../../inst
	cd "$(BUILD_DIR)" && cp nvbio-test/nvbio-test ../../../inst

clean:
	(cd "$(BUILD_DIR)" && rm -f CMakeCache.txt && ($(MAKE) clean -f Makefile))
	cd "$(BUILD_DIR)" && rm -f *.so
	cd "$(BUILD_DIR)" && rm -f *.o
	cd "$(BUILD_DIR)" && rm -f nvBWT nvBowtie nvbio-test
	cd ../../../inst && rm -f nvBWT nvBowtie nvbio-test
