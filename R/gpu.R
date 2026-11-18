#' Is a usable NVIDIA GPU available on this machine?
#'
#' Internal helper used by the vignette and by the unit tests. Building the
#' package does not need a GPU, but running the aligner does, so any R code
#' that calls nvBowtie, nvBWT or nvbio_tests has to be able to tell whether it
#' is on a build-only host (CI runner, r-universe, cluster login node,
#' container) or on a machine where the examples can really run.
#'
#' Two levels are checked, and both are validated rather than trusted:
#' \code{nvidia-smi} prints its own failures ("Failed to initialize NVML: Unknown
#' Error", "NVIDIA-SMI has failed because it couldn't communicate with the
#' NVIDIA driver") as ordinary output, sometimes on stdout, so the return code
#' and the text of the answer are both examined.
#'
#' @param min_free_mb Optional. When given, a GPU only counts as usable if at
#'   least one device reports that much free device memory. This covers the
#'   case of a shared GPU builder where the card is present but already fully
#'   allocated by someone else: nvbio then dies at context creation with
#'   "cuda-memory-check: out of memory", which is not a bug in this package.
#' @return `TRUE` or `FALSE`. Never raises an error, so it is safe to call
#'   unconditionally.
#' @author Franck RICHARD
#' @keywords internal
has_gpu <- function(min_free_mb = NULL) {
    smi <- Sys.which("nvidia-smi")
    if (!nzchar(smi) && file.exists("/usr/bin/nvidia-smi")) {
        smi <- "/usr/bin/nvidia-smi"
    }
    if (!nzchar(smi) || !file.exists(smi)) {
        return(FALSE)
    }

    # Device visibility: "GPU 0: NVIDIA ..." one line per device.
    listed <- tryCatch(
        suppressWarnings(system2(smi, args = "-L", stdout = TRUE, stderr = FALSE)),
        error = function(e) character()
    )
    if (!is.null(attr(listed, "status"))) {
        return(FALSE)
    }
    if (!any(grepl("^GPU [0-9]+:", listed))) {
        return(FALSE)
    }

    if (is.null(min_free_mb)) {
        return(TRUE)
    }

    # Free device memory, largest across devices. nounits + no header keeps the
    # answer a bare number; anything else is diagnostic text, so it is filtered.
    free <- tryCatch(
        suppressWarnings(system2(
            smi,
            args = "--query-gpu=memory.free --format=csv,noheader,nounits",
            stdout = TRUE, stderr = FALSE
        )),
        error = function(e) character()
    )
    if (!is.null(attr(free, "status"))) {
        return(FALSE)
    }
    values <- suppressWarnings(as.numeric(grep("^[0-9]+$", free, value = TRUE)))
    if (!length(values)) {
        return(FALSE)
    }
    max(values) >= min_free_mb
}

#' How much free device memory does a given entry point need?
#'
#' Internal helper backing the \code{\\examplesIf} guards and the unit test
#' skip. Building never needs a GPU, but a CUDA context alone takes a few
#' hundred MB of device memory, and Bioconductor's GPU builders share their card
#' with other tenants: amarone reported 45369 of 46068 MiB already in use, which
#' is why nvBWT died at context creation there while biocgpu, with an idle card,
#' was fine. The floors below are what those entry points were seen to need,
#' with margin; they are not measurements of a peak. They can be overridden
#' without reinstalling, with \code{options(RbowtieCuda.minFreeMemMB = 8192)}.
#'
#' @param which \code{"light"} for nvBWT()/nvBowtie() on the bundled lambda
#'   phage data, \code{"tests"} for nvbio_tests(), which walks cudaMalloc up to
#'   1 GB and then runs the alignment test on top of a context.
#' @return A single number of megabytes.
#' @author Franck RICHARD
#' @keywords internal
gpu_min_free_mb <- function(which = c("light", "tests")) {
    override <- getOption("RbowtieCuda.minFreeMemMB")
    if (is.numeric(override) && length(override) == 1L &&
        is.finite(override) && override >= 0) {
        return(as.numeric(override))
    }
    switch(match.arg(which), light = 2048, tests = 4096)
}

#' Gate used by the Rd examples of nvBWT(), nvBowtie() and nvbio_tests()
#'
#' Internal helper, called from the \code{\\examplesIf} condition. It answers the
#' question and, when the answer is no, says so: an example that is skipped
#' leaves nothing at all in the \code{R CMD check} log, which on a shared GPU
#' builder reads like the documentation is broken rather than like the card was
#' fully allocated by someone else.
#'
#' @param which passed to \code{gpu_min_free_mb}.
#' @return \code{TRUE} when the example can run, \code{FALSE} otherwise.
#' @author Franck RICHARD
#' @keywords internal
gpu_example_gate <- function(which = c("light", "tests")) {
    which <- match.arg(which)
    need <- gpu_min_free_mb(which)
    if (has_gpu(min_free_mb = need)) {
        return(TRUE)
    }
    free <- gpu_free_mb()
    what <- switch(which, light = "nvBWT()/nvBowtie()", tests = "nvbio_tests()")
    message(paste0(
        "Skipping the ", what, " example: it needs about ", need,
        " MB of free device memory, and ",
        if (is.finite(free)) {
            paste0("the largest device here reports ",
                   format(free, big.mark = " "), " MB free.")
        } else {
            "no usable GPU was found on this machine."
        }
    ))
    FALSE
}

#' Free device memory of the largest device, in MB
#'
#' Internal helper, used only to say \emph{why} an example or a test was
#' skipped. Returns \code{NA_real_} when the number cannot be read.
#'
#' @return A number of megabytes, or \code{NA_real_}.
#' @author Franck RICHARD
#' @keywords internal
gpu_free_mb <- function() {
    smi <- Sys.which("nvidia-smi")
    if (!nzchar(smi) && file.exists("/usr/bin/nvidia-smi")) {
        smi <- "/usr/bin/nvidia-smi"
    }
    if (!nzchar(smi) || !file.exists(smi)) {
        return(NA_real_)
    }
    free <- tryCatch(
        suppressWarnings(system2(
            smi,
            args = "--query-gpu=memory.free --format=csv,noheader,nounits",
            stdout = TRUE, stderr = FALSE
        )),
        error = function(e) character()
    )
    if (!is.null(attr(free, "status"))) {
        return(NA_real_)
    }
    values <- suppressWarnings(as.numeric(grep("^[0-9]+$", free, value = TRUE)))
    if (!length(values)) {
        return(NA_real_)
    }
    max(values)
}
