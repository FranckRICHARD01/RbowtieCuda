test_alignement <- function() {
    # nvbio_tests() runs the GPU test binary. A build/check machine is not
    # required to have a GPU (r-universe, CI runners, cluster login nodes), so
    # skip the test rather than fail when no device can be reached.
    # Bioconductor's GPU builders (biocgpu, kakapo1, amarone) still run it.
    if (!RbowtieCuda:::has_gpu()) {
        message("No usable GPU available: skipping nvbio_tests()")
        return(invisible(TRUE))
    }
    checkEquals(nvbio_tests(), "nvbio test passed.")
}
