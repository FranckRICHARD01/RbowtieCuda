nvbio_module(nvBWT)

find_package(CUDAToolkit REQUIRED)
find_package(OpenMP REQUIRED)

enable_language(CUDA)

addsources(
nvBWT.cu
filelist.cpp
)

add_library(nvBWT STATIC ${nvBWT_srcs})
target_link_libraries(nvBWT PRIVATE nvbio zlibstatic crcstatic moderngpu ${SYSTEM_LINK_LIBRARIES} CUDA::cuda_driver CUDA::cudart)
set_target_properties(nvBWT PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if(OpenMP_CXX_FOUND)
    target_link_libraries(nvBWT PRIVATE OpenMP::OpenMP_CXX)
    target_compile_options(nvBWT PRIVATE ${OpenMP_CXX_FLAGS})
endif()
