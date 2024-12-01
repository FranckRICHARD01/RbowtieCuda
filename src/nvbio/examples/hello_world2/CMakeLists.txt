nvbio_module(hello_world2)

addsources(
hello_world2.cu
)

cuda_add_executable(hello_world2 ${hello_world2_srcs})
target_link_libraries(hello_world2 nvbio zlibstatic lz4 crcstatic ${SYSTEM_LINK_LIBRARIES})
