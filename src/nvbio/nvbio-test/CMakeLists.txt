nvbio_module(nvbio-test)

find_package(CUDAToolkit REQUIRED)
find_package(OpenMP REQUIRED)

enable_language(CUDA)

addsources(
alignment_test.cu
alloc_test.cu
bloom_filter_test.cu
bwt_test.cpp
cache_test.cpp
condtion_test.cu
fasta_test.cpp
fastq_test.cpp
fmindex_test.cu
nvbio-test.cpp
packedstream_test.cpp
qgram_test.cu
rank_test.cu
string_set_test.cu
sum_tree_test.cpp
syncblocks_test.cu
utils.h
work_queue_test.cu
sequence_test.cu
wavelet_test.cu
)

add_executable(nvbio-test ${nvbio-test_srcs})
target_link_libraries(nvbio-test PRIVATE nvbio zlibstatic crcstatic lz4 ${SYSTEM_LINK_LIBRARIES} CUDA::cuda_driver CUDA::cudart)
set_target_properties(nvbio-test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if(OpenMP_CXX_FOUND)
    target_link_libraries(nvbio-test PRIVATE OpenMP::OpenMP_CXX)
    target_compile_options(nvbio-test PRIVATE ${OpenMP_CXX_FLAGS})
endif()

