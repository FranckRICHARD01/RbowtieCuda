nvbio_module(hello_world)

addsources(
hello_world.cu
)

cuda_add_executable(hello_world ${hello_world_srcs})
target_link_libraries(hello_world nvbio zlibstatic lz4 crcstatic ${SYSTEM_LINK_LIBRARIES})
