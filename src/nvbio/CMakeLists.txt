cmake_minimum_required(VERSION 3.22)

message(WARNING "CUDACXX=$ENV{CUDACXX}")
message(WARNING "PATH=$ENV{PATH}")

if (UNIX AND NOT APPLE)
    # Retrieve glibc version using ldd
    execute_process(
        COMMAND ldd --version
        OUTPUT_VARIABLE LDD_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Extract major and minor version numbers
    string(REGEX MATCH "ldd \\(.*\\) ([0-9]+)\\.([0-9]+)" _ UNUSED "${LDD_OUTPUT}")
    set(GLIBC_MAJOR ${CMAKE_MATCH_1})
    set(GLIBC_MINOR ${CMAKE_MATCH_2})

    # Display a warning if glibc >= 2.41
    if(GLIBC_MAJOR GREATER 2 OR (GLIBC_MAJOR EQUAL 2 AND GLIBC_MINOR GREATER_EQUAL 41))
        message(WARNING "
Attention, the glibc library on your distribution is greater than or equal to 2.41.
An incompatibility with the Cuda library currently prevents compilation of this package with this version of glibc.
If you still want to address this issue before NVIDIA releases a patch, you can follow this method:
https://github.com/FranckRICHARD01/RbowtieCuda/issues/2
        ")
    endif()
endif()

### IMPROVED COMPILER DETECTION BLOCK ###

if(NOT DEFINED CMAKE_CUDA_COMPILER)
    # Priority 1: Explicit environment variable
    if(DEFINED ENV{CUDACXX})
        set(CMAKE_CUDA_COMPILER "$ENV{CUDACXX}" CACHE FILEPATH "CUDA compiler")
    else()
        # Priority 2: Automatic search for known compilers
        message(WARNING "CUDACXX not set. Searching for 'nvcc' or 'nvc++'...")
        find_program(CUDA_COMPILER_PATH
            NAMES nvcc nvc++   # <-- KEY CHANGE: Search for both
            HINTS
                # Common paths for nvcc
                /usr/local/cuda/bin
                /opt/nvidia/hpc_sdk/Linux_x86_64/*/cuda/*/bin
                # Common paths for nvc++ (HPC SDK)
                /opt/nvidia/hpc_sdk/Linux_x86_64/*/compilers/bin
                # Generic path for installations via modules
                /software/*/bin
        )

        if(CUDA_COMPILER_PATH)
            set(CMAKE_CUDA_COMPILER "${CUDA_COMPILER_PATH}" CACHE FILEPATH "CUDA compiler")
            message(WARNING "Found CUDA compiler: ${CMAKE_CUDA_COMPILER}")
        else()
            # Error if neither is found
            message(FATAL_ERROR "CUDA compiler (nvcc or nvc++) not found. Please set the CUDACXX environment variable or add the compiler's directory to your PATH.")
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# CUDA host compiler probe
#
# nvcc accepts only a bounded range of host (g++) versions: a recent distro
# (e.g. Ubuntu 26.04 with GCC 15) combined with a CUDA 12.4 toolkit makes
# CMake's own CUDA compiler identification fail before any of our code runs.
# Instead of guessing version tables, compile one tiny .cu with nvcc; if the
# default host compiler is rejected, use the newest installed g++-NN that nvcc
# accepts. This runs BEFORE project() because that is where identification
# happens, and it never overrides a choice made by the user or the builder.
# ---------------------------------------------------------------------------
if(UNIX AND NOT APPLE AND NOT DEFINED CMAKE_CUDA_HOST_COMPILER AND NOT DEFINED ENV{CMAKE_CUDA_HOST_COMPILER}
   AND NOT CMAKE_CUDA_COMPILER MATCHES "nvc\\+\\+$")
    set(_nvbio_probe_dir "${CMAKE_CURRENT_BINARY_DIR}/_nvbio_host_probe")
    file(MAKE_DIRECTORY "${_nvbio_probe_dir}")
    file(WRITE "${_nvbio_probe_dir}/host_probe.cu"
         "#include <cuda_runtime.h>\n__global__ void nvbio_host_probe() {}\n")

    # nvcc's own default first, then every g++-NN installed on the host, newest first.
    set(_nvbio_host_candidates "nvcc-default")
    foreach(_nvbio_gcc_ver RANGE 22 9 -1)
        find_program(_nvbio_gpp_candidate NAMES "g++-${_nvbio_gcc_ver}")
        if(_nvbio_gpp_candidate)
            list(APPEND _nvbio_host_candidates "${_nvbio_gpp_candidate}")
        endif()
    endforeach()

    set(_nvbio_host_ok "")
    set(_nvbio_probe_log "")
    foreach(_nvbio_host IN LISTS _nvbio_host_candidates)
        if(_nvbio_host STREQUAL "nvcc-default")
            execute_process(
                COMMAND "${CMAKE_CUDA_COMPILER}" -c "${_nvbio_probe_dir}/host_probe.cu" -o "${_nvbio_probe_dir}/host_probe.o"
                RESULT_VARIABLE _nvbio_probe_rc OUTPUT_VARIABLE _nvbio_probe_out ERROR_VARIABLE _nvbio_probe_err
                TIMEOUT 900)
        else()
            execute_process(
                COMMAND "${CMAKE_CUDA_COMPILER}" -ccbin "${_nvbio_host}" -c "${_nvbio_probe_dir}/host_probe.cu" -o "${_nvbio_probe_dir}/host_probe.o"
                RESULT_VARIABLE _nvbio_probe_rc OUTPUT_VARIABLE _nvbio_probe_out ERROR_VARIABLE _nvbio_probe_err
                TIMEOUT 900)
        endif()
        string(STRIP "${_nvbio_probe_err}" _nvbio_probe_err)
        string(APPEND _nvbio_probe_log "\n  host compiler '${_nvbio_host}' -> exit ${_nvbio_probe_rc}${_nvbio_probe_err}")
        if(_nvbio_probe_rc EQUAL 0)
            set(_nvbio_host_ok "${_nvbio_host}")
            break()
        endif()
    endforeach()

    if(_nvbio_host_ok AND NOT _nvbio_host_ok STREQUAL "nvcc-default")
        set(CMAKE_CUDA_HOST_COMPILER "${_nvbio_host_ok}" CACHE FILEPATH "Host compiler for nvcc, chosen by the nvbio configure probe" FORCE)
        message(WARNING
          "nvcc rejected the default host compiler; compiling device code with ${_nvbio_host_ok} "
          "(CMAKE_CUDA_HOST_COMPILER). The C++ compiler picked for the host code may differ from it: "
          "if a link or ABI mismatch appears, configure with -DCMAKE_CXX_COMPILER=${_nvbio_host_ok} as well.")
    elseif(NOT _nvbio_host_ok)
        message(WARNING
          "nvcc failed to compile a trivial .cu with every host compiler tried "
          "(nvcc: ${CMAKE_CUDA_COMPILER}).${_nvbio_probe_log}\n"
          "Keeping nvcc's default; a CUDA compiler identification error below is the expected consequence.")
    endif()
    file(REMOVE_RECURSE "${_nvbio_probe_dir}")
endif()

project(nvBio LANGUAGES CXX CUDA)

execute_process(
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../inst/scripts/check_requirements.sh
    OUTPUT_VARIABLE DIAGNOSTIC_OUTPUT
    ERROR_VARIABLE DIAGNOSTIC_ERROR
    RESULT_VARIABLE DIAGNOSTIC_RESULT
)

message(WARNING "--- CUDA Environment Diagnostic ---\n${DIAGNOSTIC_OUTPUT}\n")

if (NOT DIAGNOSTIC_RESULT EQUAL 0)
    message(FATAL_ERROR "CUDA environment check failed (code: ${DIAGNOSTIC_RESULT}).")
endif()

# *** FIX: Set CMAKE_CUDA_COMPILER unconditionally and BEFORE find_package(CUDA) ***
#set(CMAKE_CUDA_COMPILER /software/u24/nvhpc/24.7/Linux_x86_64/24.7/compilers/bin/nvcc)
#set(CMAKE_CUDA_COMPILER /opt/nvidia/hpc_sdk/Linux_x86_64/25.1/cuda/12.6/bin/nvcc)
#set(CMAKE_CUDA_COMPILER /usr/local/cuda-12.9/bin/nvcc)

# Automatic GPU architecture detection.
# nvidia-smi also writes its complaints to stdout ("NVIDIA-SMI has failed
# because it couldn't communicate with the NVIDIA driver", "Failed to
# initialize NVML: Unknown Error"), so its answer must be validated line by
# line. On a GPU-less build host that text used to flow straight into
# CMAKE_CUDA_ARCHITECTURES ("Unknown CUDA architecture specifier ..."), which
# is the second half of issue #4: fixing only check_requirements.sh is not
# enough to configure on a build-only machine.
execute_process(
  COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
  OUTPUT_VARIABLE NVBIO_SMICAP_OUTPUT
  RESULT_VARIABLE NVBIO_SMICAP_RESULT
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# A genuine answer is one "N.N" per GPU; anything else is diagnostic text.
set(GPU_ARCHS "")
if(NVBIO_SMICAP_RESULT EQUAL 0 AND NVBIO_SMICAP_OUTPUT)
    string(REPLACE "\n" ";" NVBIO_SMICAP_LINES "${NVBIO_SMICAP_OUTPUT}")
    foreach(_a IN LISTS NVBIO_SMICAP_LINES)
        string(STRIP "${_a}" _a)
        if(_a MATCHES "^[0-9]+\\.[0-9]+$")
            list(APPEND GPU_ARCHS "${_a}")
        endif()
    endforeach()
    if(GPU_ARCHS)
        list(REMOVE_DUPLICATES GPU_ARCHS)
    endif()
endif()

# Oldest CUDA target supported by the installed toolkit (drives the PTX JIT fallback).
# CUDA 13 removed Pascal/Volta/Maxwell; Turing sm_75 is the oldest valid target.
if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 13.0)
    set(NVBIO_MIN_ARCH "75")
else()
    set(NVBIO_MIN_ARCH "70")
endif()

# Native arch(es) detected on this host (may be several GPUs), normalized to NN form.
set(NVBIO_NATIVE_ARCHS "")
foreach(_a IN LISTS GPU_ARCHS)
    string(REPLACE "." "" _num "${_a}")
    list(APPEND NVBIO_NATIVE_ARCHS "${_num}")
endforeach()

if(NVBIO_NATIVE_ARCHS)
    # Native SASS for detected GPU(s) + PTX JIT fallback so the fatbin runs on any arch >= floor.
    # Skip the extra PTX when a native arch already equals the floor (its own compute_XX is embedded).
    list(FIND NVBIO_NATIVE_ARCHS "${NVBIO_MIN_ARCH}" _min_present)
    if(_min_present LESS 0)
        set(CMAKE_CUDA_ARCHITECTURES "${NVBIO_NATIVE_ARCHS}" "${NVBIO_MIN_ARCH}-virtual")
        message(WARNING "CUDA arch: native=${NVBIO_NATIVE_ARCHS} + PTX JIT fallback sm_${NVBIO_MIN_ARCH}")
    else()
        set(CMAKE_CUDA_ARCHITECTURES "${NVBIO_NATIVE_ARCHS}")
        message(WARNING "CUDA arch: native=${NVBIO_NATIVE_ARCHS} (floor already embedded)")
    endif()
else()
    # No detectable GPU (nvidia-smi unavailable/failing): oldest supported SASS + PTX fallback.
    set(CMAKE_CUDA_ARCHITECTURES "${NVBIO_MIN_ARCH}" "${NVBIO_MIN_ARCH}-virtual")
    message(WARNING "Could not detect GPU architecture. Using sm_${NVBIO_MIN_ARCH} SASS + PTX JIT fallback.")
endif()

#set(CMAKE_CXX_COMPILER_ID NVHPC)

if (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(USING_NVHPC TRUE)
else()
    set(USING_NVHPC FALSE)
endif()

#if (USING_NVHPC)
#    #set(CMAKE_C_COMPILER nvc)
#    #set(CMAKE_CXX_COMPILER nvc++)
#    set(CMAKE_CUDA_COMPILER nvc++)
#else()
#    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
#endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# warning 
# If find_package(CUDA required) has trouble finding the correct path to the cuda toolkit, uncheck this parameter and specify the correct path.
#set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda-12.9")



# Disable CUB warning
add_definitions(-DTHRUST_IGNORE_CUB_VERSION_CHECK=1)

# available options
option(CUDA_VERBOSE_PTXAS
  "Enable verbose output from PTXAS"
  ON)

option(PROFILING
  "Enable profiling"
  OFF)

option(WERROR
  "Treat compiler warnings as errors"
  OFF)

set(NVBIO_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# import our own nvbio_* macros
include("${NVBIO_SOURCE_ROOT}/cmake-local/nvbio.cmake")
# grab gcc flags
include("${NVBIO_SOURCE_ROOT}/cmake-local/gcc.cmake")

# default to release builds
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Building...

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DTHRUST_DEBUG)
  add_definitions(-DCUDA_DEBUG)
  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  #set(CUDA_NVCC_DEBUG_FLAGS "-O0 -g -lineinfo ")
  set(CUDA_NVCC_DEBUG_FLAGS "-O0 -g -G --ptxas-options=-O0")
  if (UNIX AND NOT APPLE AND NOT USING_NVHPC)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
  endif()
else()
  set(CUDA_NVCC_DEBUG_FLAGS "--ptxas-options=-O3 -O3")
  if (UNIX AND NOT APPLE AND NOT USING_NVHPC)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-old-style-definition -O3")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-enum-enum-conversion -Wno-deprecated-declarations -Wno-unused-result -O3")
  endif()
endif()


if (CMAKE_SYSTEM_PROCESSOR MATCHES "Intel" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
add_definitions(-DPLATFORM_X86)
endif()

find_package(CUDAToolkit REQUIRED)
find_package(Thrust REQUIRED)
find_package(CUB REQUIRED)
find_package(Doxygen)

#enable_language(CUDA)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set up required system libs
set(SYSTEM_LINK_LIBRARIES "")
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
  # shm_* live in librt on Linux
  list(APPEND SYSTEM_LINK_LIBRARIES rt)
endif()

if (MSVC_IDE)
  # suppress automatic regeneration of VS project files
  #   set(CMAKE_SUPPRESS_REGENERATION ON)
  # compilation options
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MD /D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
endif() 

if (APPLE AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  find_program(CLANG_BIN clang)
  set(CUDA_NVCC_FLAGS " -ccbin ${CLANG_BIN} -Xcompiler -stdlib=libstdc++")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -stdlib=libstdc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
endif()

# Optimal configuration for NVHPC
if(USING_NVHPC)
  set(BASE_NVHPC_ARCH_FLAGS "-gpu=nordc,cc70,cc75,cc80,cc86,cc89")
  
  if(GPU_ARCHS)
    # Extract specific architectures
    list(REMOVE_DUPLICATES GPU_ARCHS)
    string(REPLACE "." "" GPU_ARCHS_CLEAN ${GPU_ARCHS})
    list(TRANSFORM GPU_ARCHS_CLEAN PREPEND "cc")
    
    # Combine with base flags
    string(JOIN "," ARCH_STRING ${GPU_ARCHS_CLEAN})
    set(NVHPC_ARCH_FLAGS "${BASE_NVHPC_ARCH_FLAGS},${ARCH_STRING}")
  else()
    set(NVHPC_ARCH_FLAGS ${BASE_NVHPC_ARCH_FLAGS})
  endif() 

endif()

# Configuration for nvcc
# NOTE: superseded by CMAKE_CUDA_ARCHITECTURES above. CUDA_NVCC_FLAGS is not consumed by
# modern CMake; kept only for reference. Do NOT reintroduce sm_70 (removed in CUDA 13).
if(NOT USING_NVHPC)
  if(GPU_ARCHS)
    # Processing for nvcc
    list(REMOVE_DUPLICATES GPU_ARCHS)
    string(REPLACE "." "" GPU_ARCHS_CLEAN ${GPU_ARCHS})
    list(TRANSFORM GPU_ARCHS_CLEAN PREPEND "sm_")
    
    foreach(ARCH IN LISTS GPU_ARCHS_CLEAN)
      string(SUBSTRING ${ARCH} 3 -1 ARCH_NUM)
      list(APPEND NVCC_ARCH_FLAGS 
           -gencode arch=compute_${ARCH_NUM},code=compute_${ARCH_NUM}
           -gencode arch=compute_${ARCH_NUM},code=${ARCH})
    endforeach()
  else()
    set(NVCC_ARCH_FLAGS
      "-gencode arch=compute_70,code=sm_70
      -gencode arch=compute_75,code=sm_75
      -gencode arch=compute_80,code=sm_80
      -gencode arch=compute_86,code=sm_86
      -gencode arch=compute_89,code=sm_89")
  endif()
endif()

if (USING_NVHPC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fast")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_DEBUG_FLAGS}
    ${NVHPC_ARCH_FLAGS}
    -fastmath")
else()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_NVCC_DEBUG_FLAGS} 
    ${NVCC_ARCH_FLAGS}
    -use_fast_math")
endif()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-unused-result -Wno-deprecated-enum-enum-conversion -Wno-deprecated-declarations -Wno-deprecated-gpu-targets -diag-suppress 20011 -diag-suppress 20014 -diag-suppress 68 -diag-suppress 1650 -diag-suppress 177 -diag-suppress 186")


set(cudpp_INCLUDE_DIRS ${NVBIO_SOURCE_ROOT})

set(INCLUDE_DIRS
        ${NVBIO_SOURCE_ROOT}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NVBIO_SOURCE_ROOT}/contrib
        ${NVBIO_SOURCE_ROOT}/contrib/crc
        ${NVBIO_SOURCE_ROOT}/contrib/bamtools
        ${NVBIO_SOURCE_ROOT}/contrib/moderngpu/include
        )

include_directories(${NVBIO_SOURCE_ROOT})
include_directories("${NVBIO_SOURCE_ROOT}/contrib")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/crc")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/bamtools")
include_directories("${CUDA_TOOLKIT_ROOT_DIR}/targets/x86_64-linux/include")
include_directories(SYSTEM "${NVBIO_SOURCE_ROOT}/contrib/moderngpu/include")

# export our include directories in the cache for external projects to pick up
set(nvBio_INCLUDE_DIRECTORIES
        ${INCLUDE_DIRS} CACHE STRING "nvbio include paths")
set(nvBio_INCLUDE
        "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(contrib/zlib-1.2.7)
add_subdirectory(contrib/lz4)
add_subdirectory(contrib/crc)
add_subdirectory(contrib/bamtools)
add_subdirectory(contrib/moderngpu)
#add_subdirectory(contrib/htslib)

add_subdirectory(nvbio)
add_subdirectory(nvBowtie) 
#add_subdirectory(nvFM-server)
add_subdirectory(nvBWT)
#add_subdirectory(nvSetBWT)
#add_subdirectory(nvSSA)
#add_subdirectory(nvExtractReads)
add_subdirectory(nvbio-test) 
#add_subdirectory(nvbio-aln-diff)
#add_subdirectory(nvMicroAssembly)
#add_subdirectory(sufsort-test)
#add_subdirectory(sw-benchmark)

##add_subdirectory(examples/waveletfm)
#add_subdirectory(examples/proteinsw)
#add_subdirectory(examples/seeding)
#add_subdirectory(examples/fmmap)
##add_subdirectory(examples/qmap)
#add_subdirectory(examples/mem)

if(USING_NVHPC)
  set(ARCH_FLAGS "Detected NVHPC Compiler. NVHPC Flags:        ${NVHPC_ARCH_FLAGS}")
else()
  set(ARCH_FLAGS "Detected NVCC Compiler. NVCC Flags:         ${NVCC_ARCH_FLAGS}")
endif()

message(WARNING "================================================\n"
                "CUDA Configuration Summary\n"
                "Compiler:            ${CMAKE_CUDA_COMPILER}\n"
                "Build type:          ${CMAKE_BUILD_TYPE}\n"
                "GPU Architectures:   ${GPU_ARCHS}\n"
                "${ARCH_FLAGS}\n"
                "================================================")
#Install targets
install(TARGETS ${NVBIO_BIN_LIST} DESTINATION bin)
install(FILES ${NVBIO_WRAPPER_SCRIPTS}
    PERMISSIONS OWNER_REAK OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

#nvbio_module(nvbio-project)
nvbio_doxygen()