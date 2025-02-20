nvbio_module(nvBowtie)

find_package(CUDAToolkit REQUIRED)
find_package(OpenMP REQUIRED)

enable_language(CUDA)

addsources(
nvBowtie.cu
)

nvbio_add_module_directory(bowtie2)
nvbio_add_module_directory(bowtie2/cuda)

add_executable(nvBowtie ${nvBowtie_srcs})
target_link_libraries(nvBowtie PRIVATE nvbio zlibstatic lz4 crcstatic ${SYSTEM_LINK_LIBRARIES} CUDA::cuda_driver CUDA::cudart)
set_target_properties(nvBowtie PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if(OpenMP_CXX_FOUND)
    target_link_libraries(nvBowtie PRIVATE OpenMP::OpenMP_CXX)
    target_compile_options(nvBowtie PRIVATE ${OpenMP_CXX_FLAGS})
endif()
