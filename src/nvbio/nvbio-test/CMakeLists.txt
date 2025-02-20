nvbio_module(nvbio-test)

find_package(CUDAToolkit REQUIRED)
find_package(OpenMP REQUIRED)

enable_language(CUDA)

addsources(
alignment_test.cu
alloc_test.cu
bloom_filter_test.cu
bwt_test.cu
cache_test.cu
condtion_test.cu
fasta_test.cu
fastq_test.cu
fmindex_test.cu
nvbio-test.cu
packedstream_test.cu
qgram_test.cu
rank_test.cu
string_set_test.cu
sum_tree_test.cu
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

