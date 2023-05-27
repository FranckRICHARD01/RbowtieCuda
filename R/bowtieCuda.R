globalVariables(".")

#' @name nvBWT
#' @title nvBWT : builds the BWT indices of the reference FASTA files
#' @description This function can be use to call wrapped \code{nvBWT}
#' binary.
#' @param myinput \code{Character} vector. A path to a .fa file
#' @param output \code{Character} vector. A path to a index file
#' used for the alignment output.
#' @param options \code{Character} vector. Options
#' Example: " --verbosity  "
#'
#' @details
#' Available options:
#'
#'    -v       | --verbosity     int (0-6) [5]     // select the verbosity level
#'    -m       | --max-length    int       [inf]   // clamp input length
#'    -b       | --byte-packing  [default]         // output a byte-encoded .pac file
#'    -w       | --word-packing                    // output a word-encoded .wpac file (more efficient)
#'    -c       | --crc                             // compute CRCs
#'    -d       | --device                          // select a cuda device
#'
#' @author Franck RICHARD
#' @return An invisible \code{Integer} of call
#' status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export nvBWT
#' @examples
#' td <- tempdir()
#' fa_file <- system.file(package="RbowtieCuda", "extdata", "bt2", "refs", "lambda_virus.fa")
#' nvBWT(myinput=fa_file, output=file.path(td, "index"), options="")


nvBWT <- function(myinput, output, options = NULL) {
    if (R.Version()$arch == "i386") {
        return("nvBWT is not available for 32bit, please use 64bit R instead")
    }

    myinput <- file.path(tools::file_path_as_absolute(dirname(myinput)), basename(myinput))
    myinput <- trimws(as.character(myinput))
    myinput <- shQuote(myinput)
    #print(paste("myinput: ", myinput))

    output <- file.path(tools::file_path_as_absolute(dirname(output)), basename(output))
    output <- trimws(as.character(output))
    output <- shQuote(output)
    #print(paste("output: ", output))

    if (is.null(options))
        args1 <- paste(myinput, " ", output)
    else
        args1 <- paste(myinput, " ", output, options)

    .callbinary2(
        bin1 = "nvBWT",
        args1
    )
}


#' @name nvBowtie
#' @title nvBowtie align function
#' @description nvBowtie is a GPU-accelerated re-engineering of Bowtie2, a very widely used short-read aligner.
#' @param index \code{Character} vector. Index file created by nvBWT.
#' @param output_file \code{Character} vector. The alignment output file (.sam o .bam)
#' @param options \code{Character} vector. Indicate additional options here to modify the processing. The complete
#'  list of available options can be accessed by typing \code{nvBowtie_usage()}.
#' @param seq1 \code{Character} vector. For single-end sequencing, it contains sequence file paths.
#' For paired-end sequencing, it can be file path with #1 mates paired with file paths in seq2.
#' @param seq2 \code{Character} vector. It contains file paths with #2 mates paired with file path in seq1.
#' For single-end sequencoing files, it must be \code{NULL}.'
#' @return An invisible \code{Integer} of call
#' status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export nvBowtie
#' @examples
#' td <- tempdir()
#'
#' ## Building index
#' fa_file <- system.file(package="RbowtieCuda", "extdata", "bt2", "refs", "lambda_virus.fa")
#' nvBWT(myinput=fa_file, output=file.path(td, "index"), options="")
#'
#' ## Alignments
#' read_1 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_1.fastq")
#' read_2 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_2.fastq")
#'
#' ## Sam file created with paired-end files (paired mates are here reverse-forward)
#' nvBowtie(file.path(td, "index"), file.path(td, "my_result.sam"), options="--rf", seq1=read_1, seq2=read_2)
#'
#' ## Bam file created with single-end file
#' nvBowtie(file.path(td, "index"), file.path(td, "my_result.bam"), options="", seq1=read_1)



nvBowtie <- function(index, output_file, options, seq1, seq2 = NULL) {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit, 
        please use 64bit R instead")
    }

    index <- file.path(tools::file_path_as_absolute(dirname(index)), basename(index))
    index <- trimws(as.character(index))
    index <- shQuote(index)
    #print(paste("index: ", index))

    output_file <- file.path(tools::file_path_as_absolute(dirname(output_file)), basename(output_file))
    output_file <- trimws(as.character(output_file))
    output_file <- shQuote(output_file)
    
    #print(paste("output_file: ", output_file))


    if (!is.null(seq1))
      seq1 <- trimws(as.character(seq1))

    # If paired mates are provided then they should be the same length
    if (!is.null(seq2)) {
      seq2 <- trimws(as.character(seq2))
      if (length(seq1) != length(seq2)) {
        stop("The lengths of arguments ",
             "`seq1` and `seq2` should be the same length")
      }
    }

    if (!is.null(seq2)) {
      args1 <- paste(options, "--file-ref -x ", index, " -1 ", seq1,
      " -2 ", seq2, " -S ", output_file)
    } else {
      args1 <- paste(options, "--file-ref -x ", index, " -U ", seq1,
      " -S ", output_file)
    }

    #print(paste("nvBowtie ", args1))

    .callbinary2(
        bin1 = "nvBowtie",
        args1
    )
}



#' @name nvBowtie_version
#' @title Print version information
#' @description Calling nvBowtie_version() prints the version information of 
#' the RBowtieCuda package used.
#' @author Franck RICHARD
#' @return An invisible \code{Integer} of call status.
#' The value is 0 when there is not selected
#' @export nvBowtie_version
#' @examples
#' nvBowtie_version()

nvBowtie_version <- function() {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit,
         please use 64bit R instead")
    }

    .callbinary2(
        bin1 = "nvBowtie",
        args1 = "--version"
    )
}

#' @name nvBowtie_usage
#' @title Print available arguments that can be passed to nvBowtie()
#' @description Calling nvBowtie_usage() prints the available arguments that can
#' be passed to the ... argument 'options' of the nvBowtie() function of the package.
#' @author Franck RICHARD
#' @return Information about available arguments that can be passed to nvBowtie().
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export nvBowtie_usage
#' @examples
#' nvBowtie_usage()



nvBowtie_usage <- function() {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit,
         please use 64bit R instead")
    }

    #invisible(.callbinary2(
    #    bin1 = "nvBowtie",
    #args1 = "-h"))

    print("options:")
    print("General:")
    print("  --verbosity         int [5]        verbosity level")
    print("  --upto       | -u   int [-1]       maximum number of reads to process")
    print("  --trim3      | -3   int [0]        trim the first N bases of 3'")
    print("  --trim5      | -5   int [0]        trim the first N bases of 5'")
    print("  --nofw                             do not align the forward strand")
    print("  --norc                             do not align the reverse-complemented strand")
    print("  --device            int [0]        select the given cuda device(s) (e.g. --device 0 --device 1 ...)")
    print("  --file-ref                         load reference from file")
    print("  --server-ref                       load reference from server")
    print("  --phred33                          qualities are ASCII characters equal to Phred quality + 33")
    print("  --phred64                          qualities are ASCII characters equal to Phred quality + 64")
    print("  --solexa-quals                     qualities are in the Solexa format")
    print("  --rg-id             string         add the RG-ID field of the SAM output header")
    print("  --rg                string,val     add an RG-TAG field of the SAM output header")
    print("Paired-End:")
    print("  --ff                               paired mates are forward-forward")
    print("  --fr                               paired mates are forward-reverse")
    print("  --rf                               paired mates are reverse-forward")
    print("  --rr                               paired mates are reverse-reverse")
    print("  --minins            int [0]        minimum insert length")
    print("  --maxins            int [500]      maximum insert length")
    print("  --overlap                          allow overlapping mates")
    print("  --dovetail                         allow dovetailing mates")
    print("  --no-mixed                         only report paired alignments")
    print("  --ungapped-mates | -ug             perform ungapped mate alignment")
    print("Seeding:")
    print("  --seed-len   | -L   int   [22]     seed lengths")
    print("  --seed-freq  | -i   {G|L|S},x,y    seed interval, as x + y*func(read-len) (G=log,L=linear,S=sqrt)")
    print("  --max-hits          int   [100]    maximum amount of seed hits")
    print("  --max-reseed | -R   int   [2]      number of reseeding rounds")
    print("Extension:")
    print("  --all        | -a                  perform all-mapping (i.e. find and report all alignments)")
    print("  --local                            perform local alignment")
    print("  --rand                             randomized seed selection")
    print("  --no-rand                          do not randomize seed hit selection")
    print("  --max-dist          int [15]       maximum edit distance")
    print("  --max-effort-init   int [15]       initial maximum number of consecutive extension failures")
    print("  --max-effort | -D   int [15]       maximum number of consecutive extension failures")
    print("  --min-ext           int [30]       minimum number of extensions per read")
    print("  --max-ext           int [400]      maximum number of extensions per read")
    print("  --fast                             apply the fast presets")
    print("  --very-fast                        apply the very-fast presets")
    print("  --sensitive                        apply the sensitive presets")
    print("  --very-sensitive                   apply the very-sensitive presets")
    print("  --fast-local                       apply the fast presets")
    print("  --very-fast-local                  apply the very-fast presets")
    print("  --sensitive-local                  apply the sensitive presets")
    print("  --very-sensitive-local             apply the very-sensitive presets")
    print("Scoring:")
    print("  --scoring           {sw|ed}        Smith-Waterman / Edit-Distance scoring")
    print("  --score-min         {G|L|S},x,y    minimum score function, as x + y*func(read-len)")
    print("  --ma                int            match bonus")
    print("  --mp                int,int        mismatch min/max penalties")
    print("  --np                int            N penalty")
    print("  --rdg               int,int        read open/extension gap penalties")
    print("  --rfg               int,int        reference open/extension gap penalties")
    print("Reporting:")
    print("  --mapQ-filter | -Q  int [0]        minimum mapQ threshold")
}



#' @name .callbinary2
#' @title Make system call for binaries
#' @description Function that makes a system call for the RbowtieCuda binaries.
#'  Note it is not intended to be used outside of the package.
#' @author Franck RICHARD
#' @param bin1 \code{Character}. The binary used for the system call.
#' @param args1 \code{Character}. The arguments to pass to the binary.
#' @param path \code{Character} Optional: If passed to function, returns
#'  the path.
#' @export .callbinary2
#' @return The output of the system call or the path provided.

.callbinary2 <- function(bin1, args1, path = NULL) {

    bin1 <- paste(file.path(system.file(package = "RbowtieCuda"), bin1))

    result <- system2(bin1, args1)

    if (!is.null(path))
        return(path)
    else
        return(result)
}