project(nvBio)

cmake_minimum_required(VERSION 3.22)
#set(CMAKE_CXX_STANDARD 03)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "-std=c++14") # something may error caused by nvcc

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

set(GPU_ARCHITECTURE "sm_60" CACHE STRING "Target GPU architecture")

set(NVBIO_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# import our own nvbio_* macros
include("${NVBIO_SOURCE_ROOT}/cmake-local/nvbio.cmake")
# grab gcc flags
include("${NVBIO_SOURCE_ROOT}/cmake-local/gcc.cmake")


# default to release builds
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Building ${CMAKE_BUILD_TYPE}")

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-DTHRUST_DEBUG)
  add_definitions(-DCUDA_DEBUG)
  set(CUDA_NVCC_DEBUG_FLAGS "-std=c++14 -g -G")
else()
  set(CUDA_NVCC_DEBUG_FLAGS "-std=c++14 -Wno-deprecated-declarations")
  if (UNIX AND NOT APPLE)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  endif()
endif()

if (CMAKE_SYSTEM_PROCESSOR MATCHES "Intel" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
add_definitions(-DPLATFORM_X86)
endif()

find_package(CUDA)
find_package(Doxygen)

find_package(OpenMP)
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
  #set(CMAKE_SUPPRESS_REGENERATION ON)
endif()

if (APPLE AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  find_program(CLANG_BIN clang)
  set(CUDA_NVCC_FLAGS " -ccbin ${CLANG_BIN} -Xcompiler -stdlib=libstdc++")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -stdlib=libstdc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
endif()

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_DEBUG_FLAGS} --gpu-architecture=${GPU_ARCHITECTURE} -use_fast_math -O3)

set(cudpp_INCLUDE_DIRS ${NVBIO_SOURCE_ROOT})

set(INCLUDE_DIRS
        ${NVBIO_SOURCE_ROOT}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NVBIO_SOURCE_ROOT}/contrib
        #${CMAKE_BINARY_DIR}/contrib/zlib-1.2.7
        ${NVBIO_SOURCE_ROOT}/contrib/crc
        ${NVBIO_SOURCE_ROOT}/contrib/bamtools
        ${NVBIO_SOURCE_ROOT}/contrib/cub
        ${NVBIO_SOURCE_ROOT}/contrib/moderngpu/include
        )

include_directories(${NVBIO_SOURCE_ROOT})
include_directories("${NVBIO_SOURCE_ROOT}/contrib")
#include_directories("${CMAKE_BINARY_DIR}/contrib/zlib-1.2.7")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/crc")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/thrust")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/bamtools")
include_directories("${NVBIO_SOURCE_ROOT}/contrib/cub")
include_directories(SYSTEM "${NVBIO_SOURCE_ROOT}/contrib/moderngpu/include")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -iquote ${CMAKE_CURRENT_SOURCE_DIR}")

# export our include directories in the cache for external projects to pick up
set(nvBio_INCLUDE_DIRECTORIES
        ${INCLUDE_DIRS} CACHE STRING "nvbio include paths")
set(nvBio_INCLUDE
        "${CMAKE_CURRENT_SOURCE_DIR}")


add_subdirectory(contrib/zlib-1.2.7)
add_subdirectory(contrib/thrust)
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
#add_subdirectory(nvbio-test)
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

install(TARGETS ${NVBIO_BIN_LIST} DESTINATION bin)
install(FILES ${NVBIO_WRAPPER_SCRIPTS}
    PERMISSIONS OWNER_REAK OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION bin)

nvbio_module(nvbio-project)
nvbio_doxygen()

