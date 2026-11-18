test_alignement <- function() {
    # nvbio_tests() runs the GPU test binary. A build/check machine is not
    # required to have a GPU (r-universe, CI runners, cluster login nodes), so
    # skip the test rather than fail when no device can be reached.
    # Bioconductor's GPU builders (biocgpu, kakapo1, amarone) still run it.
    #
    # A visibility check alone is not enough on a shared card: nvbio-test
    # cudaMallocs up to 1 GB and then runs the alignment test on top of a CUDA
    # context, which itself takes a few hundred MB. On amarone (2026-11-17)
    # 45369 of 46068 MiB of the L40S were already held by another tenant, so the
    # test started, crawled (139 ms for a 1 MB cudaMalloc) and ended in
    # "cudaErrorMemoryAllocation: out of memory". That is device contention, not
    # a regression in this package, so skip when the card cannot host the run.
    # The floor can be raised without reinstalling, with
    # options(RbowtieCuda.minFreeMemMB = 8192).
    need <- RbowtieCuda:::gpu_min_free_mb("tests")
    if (!RbowtieCuda:::has_gpu(min_free_mb = need)) {
        free <- RbowtieCuda:::gpu_free_mb()
        message(paste0(
            "Skipping nvbio_tests(): this test needs about ", need,
            " MB of free device memory, and ",
            if (is.finite(free)) {
                paste0("the largest device here reports ",
                       format(free, big.mark = " "), " MB free.")
            } else {
                "no usable GPU was found on this machine."
            }
        ))
        return(invisible(TRUE))
    }
    checkEquals(nvbio_tests(), "nvbio test passed.")
}
