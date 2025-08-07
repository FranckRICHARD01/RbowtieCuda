cmake_minimum_required(VERSION 3.22)

project(nvBio)

execute_process(
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../inst/scripts/check_requirements.sh
    OUTPUT_VARIABLE DIAGNOSTIC_OUTPUT
    ERROR_VARIABLE DIAGNOSTIC_ERROR
    RESULT_VARIABLE DIAGNOSTIC_RESULT
)

message(WARNING "--- CUDA Environment Diagnostic ---\n${DIAGNOSTIC_OUTPUT}\n")

# *** FIX: Set CMAKE_CUDA_COMPILER unconditionally and BEFORE find_package(CUDA) ***
#set(CMAKE_CUDA_COMPILER /software/u24/nvhpc/24.7/Linux_x86_64/24.7/compilers/bin/nvcc)
#set(CMAKE_CUDA_COMPILER /opt/nvidia/hpc_sdk/Linux_x86_64/25.1/cuda/12.6/bin/nvcc)
#set(CMAKE_CUDA_COMPILER /usr/local/cuda-12.8/bin/nvcc)
#set(CMAKE_CUDA_COMPILER /usr/local/cuda-12.5/bin/nvcc)

find_package(CUDAToolkit REQUIRED)

#set(CMAKE_CXX_COMPILER_ID NVHPC)

if (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(USING_NVHPC TRUE)
else()
    set(USING_NVHPC FALSE)
endif()

if (USING_NVHPC)
    #set(CMAKE_C_COMPILER nvc)
    #set(CMAKE_CXX_COMPILER nvc++)
    set(CMAKE_CUDA_COMPILER nvc++)
else()
    set(CMAKE_CUDA_COMPILER nvcc)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# warning 
# If find_package(CUDA required) has trouble finding the correct path to the cuda toolkit, uncheck this parameter and specify the correct path.
#set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda-12.8")

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
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

set(GPU_ARCHITECTURE "sm_60" CACHE STRING "Target GPU architecture")
set(COMPUTE_ARCHITECTURE "compute_60" CACHE STRING "Target COMPUTE architecture")

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
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-enum-enum-conversion -Wno-deprecated-declarations -O3")
  endif()
endif()


if (CMAKE_SYSTEM_PROCESSOR MATCHES "Intel" OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86")
add_definitions(-DPLATFORM_X86)
endif()

find_package(CUDAToolkit REQUIRED)
find_package(Thrust REQUIRED)
find_package(CUB REQUIRED)
find_package(Doxygen)

enable_language(CUDA)

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

# Automatic GPU architecture detection
execute_process(
  COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
  OUTPUT_VARIABLE GPU_ARCHS
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

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
if(NOT USING_NVHPC)
  if(GPU_ARCHS)
    # Processing for nvcc
    list(REMOVE_DUPLICATES GPU_ARCHS)
    string(REPLACE "." "" GPU_ARCHS_CLEAN ${GPU_ARCHS})
    list(TRANSFORM GPU_ARCHS_CLEAN PREPEND "sm_")
    
    foreach(ARCH IN LISTS GPU_ARCHS_CLEAN)
      string(SUBSTRING ${ARCH} 3 -1 ARCH_NUM)
      list(APPEND NVCC_ARCH_FLAGS "-gencode=arch=compute_${ARCH_NUM},code=${ARCH}")
    endforeach()
  else()
    set(NVCC_ARCH_FLAGS
      "-gencode=arch=compute_70,code=sm_70
      -gencode=arch=compute_75,code=sm_75
      -gencode=arch=compute_80,code=sm_80
      -gencode=arch=compute_86,code=sm_86
      -gencode=arch=compute_89,code=sm_89")
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

nvbio_module(nvbio-project)
nvbio_doxygen()