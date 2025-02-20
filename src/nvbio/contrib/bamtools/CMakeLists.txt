project(bamtools)

cmake_minimum_required(VERSION 3.22)

find_package(CUDAToolkit REQUIRED)

enable_language(CUDA)

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${NVBIO_SOURCE_ROOT}/nvbio)

set(BAMTOOLS_SRCS
    BamAux.h
    BamReader.cpp
    BamReader.h
    BamWriter.cpp
    BamWriter.h
    BGZF.cpp
    BGZF.h
)

add_library(bamtools STATIC ${BAMTOOLS_SRCS})
target_link_libraries(bamtools PRIVATE CUDA::cuda_driver CUDA::cudart)
set_target_properties(bamtools PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

