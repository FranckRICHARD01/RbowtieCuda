#' Is a usable NVIDIA GPU available on this machine?
#'
#' Internal helper used by the vignette and by the unit tests. Building the
#' package does not need a GPU, but running the aligner does, so any R code
#' that calls nvBowtie, nvBWT or nvbio_tests has to be able to tell whether it
#' is on a build-only host (CI runner, r-universe, cluster login node,
#' container) or on a machine where the examples can really run.
#'
#' Two levels are checked, and both are validated rather than trusted:
#' `nvidia-smi` prints its own failures ("Failed to initialize NVML: Unknown
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
