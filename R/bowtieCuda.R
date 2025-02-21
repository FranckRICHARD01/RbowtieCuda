globalVariables(".")

#' @name nvBWT
#' @title nvBWT : builds the BWT indices of the reference FASTA files
#' @description This function can be use to call wrapped \code{nvBWT}
#' binary.
#' @param myinput \code{Character} vector. A path to a .fa file
#' @param output \code{Character} vector. A path to a index file
#' used for the alignment output.
#' @param options \code{Character} vector. Options
#' Example: \code{--verbosity}
#'
#' @details
#' Available options:
#'
#' \preformatted{
#' -v       | --verbosity     int (0-6) [5]     // select the verbosity level
#' -m       | --max-length    int       [inf]   // clamp input length
#' -b       | --byte-packing  [default]         // output a byte-encoded .pac file
#' -w       | --word-packing                    // output a word-encoded .wpac file
#'                                              // (more efficient)
#' -c       | --crc                             // compute CRCs
#' -d       | --device                          // select a cuda device
#' }
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
#' fa_file <- system.file(package='RbowtieCuda', 'extdata', 'bt2', 'refs', 'lambda_virus.fa')
#' nvBWT(myinput=fa_file, output=file.path(td, 'index'), options='')


nvBWT <- function(myinput, output, options = NULL) {
    if (R.Version()$arch == "i386") {
        return("nvBWT is not available for 32bit, please use 64bit R instead")
    }

    myinput <- file.path(tools::file_path_as_absolute(dirname(myinput)),
        basename(myinput))
    myinput <- trimws(as.character(myinput))
    myinput <- shQuote(myinput)
    # show(paste('myinput: ', myinput))

    output <- file.path(tools::file_path_as_absolute(dirname(output)),
        basename(output))
    output <- trimws(as.character(output))

    unlink(paste(output, ".amb", sep = ""))
    unlink(paste(output, ".ann", sep = ""))
    unlink(paste(output, ".bwt", sep = ""))
    unlink(paste(output, ".rbwt", sep = ""))
    unlink(paste(output, ".pac", sep = ""))
    unlink(paste(output, ".rpac", sep = ""))
    unlink(paste(output, ".rsa", sep = ""))
    unlink(paste(output, ".sa", sep = ""))

    output <- shQuote(output)
    # show(paste('output: ', output))

    if (is.null(options))
        args1 <- paste(myinput, " ", output) else args1 <- paste(myinput, " ", output, options)

    .callbinary2(bin1 = "nvBWT", args1)
}

#' @name nvBowtie
#' @title nvBowtie align function
#' @description nvBowtie is a GPU-accelerated re-engineering of Bowtie2, a very widely used short-read aligner.
#' @param index \code{Character} vector. Index file created by nvBWT.
#' @param output_file \code{Character} vector. The alignment output file (.sam o .bam)
#' @param options \code{Character} vector. Specify additional options here to customize the processing. 
#' For a complete list of available options, type  \code{nvBowtie_usage()}.
#' @param seq1 \code{Character} vector. For single-end sequencing, it contains sequence file paths.
#' For paired-end sequencing, it can be file path with #1 mates paired with file paths in seq2.
#' @param seq2 \code{Character} vector. It contains file paths with #2 mates paired with file path in seq1.
#' For single-end sequencoing files, it must be \code{NULL}.
#' @return An invisible \code{Integer} of call
#' status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @details We’ve introduced several new features to nvBowtie. You can now perform alignments using the WFA
#' method by including the \code{--wfa (or --scoring wfa)} parameter. The WFA method requires a large amount 
#' of RAM on the graphics card. We therefore recommend using an Nvidia card with 8GB or more.
#' Please note that this feature is still experimental; 
#' it currently supports only end-to-end alignments and does not yet allow customization of scoring parameters. 
#' By default, it uses the following scoring: \code{match:0, mismatch:1, gap_open:1 and gap_ext:1}.
#' Additionally, the \code{--cache-writes} parameter optimizes disk write operations, resulting in faster alignments.
#' This functionality requires 4GB of RAM and is limited to paired-end alignments.
#' @references Langmead, B., & Salzberg, S. L. (2012).
#' Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export nvBowtie
#' @examples
#' td <- tempdir()
#'
#' ## Building index
#' fa_file <- system.file(package='RbowtieCuda', 'extdata', 'bt2', 'refs', 'lambda_virus.fa')
#' nvBWT(myinput=fa_file, output=file.path(td, 'index'), options='')
#'
#' ## Alignments
#' read_1 <- system.file(package='RbowtieCuda', 'extdata', 'bt2', 'reads', 'reads_1.fastq')
#' read_2 <- system.file(package='RbowtieCuda', 'extdata', 'bt2', 'reads', 'reads_2.fastq')
#'
#' ## Sam file created with paired-end files (paired mates are here reverse-forward)
#' nvBowtie(file.path(td, 'index'), file.path(td, 'my_result.sam'), options='--rf', seq1=read_1, seq2=read_2)
#'
#' ## Bam file created with single-end file
#' nvBowtie(file.path(td, 'index'), file.path(td, 'my_result.bam'), options='', seq1=read_1)



nvBowtie <- function(index, output_file, options, seq1, seq2 = NULL) {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit, please use 64bit R instead")
    }

    index <- file.path(tools::file_path_as_absolute(dirname(index)), basename(index))
    index <- trimws(as.character(index))
    index <- shQuote(index)
    # show(paste('index: ', index))

    output_file <- file.path(tools::file_path_as_absolute(dirname(output_file)),
        basename(output_file))
    output_file <- trimws(as.character(output_file))
    output_file <- shQuote(output_file)
    # show(paste('output_file: ', output_file))

    unlink(output_file)

    if (grepl("\\.fastq", options) || grepl("\\.fasta", options) || grepl("\\.fq$",
        options) || grepl("\\.fa$", options) || grepl("\\.gz$", options) || grepl("\\.sam",
        options) || grepl("\\.bam", options)) {
        stop("'Options' is a required parameter. It must always be mentioned, even if it's empty.")
    }

    if (!is.null(seq1))
        seq1 <- trimws(as.character(seq1))

    # If paired mates are provided then they should be the same length
    if (!is.null(seq2)) {
        seq2 <- trimws(as.character(seq2))
        if (length(seq1) != length(seq2)) {
            stop("The lengths of arguments ", "`seq1` and `seq2` should be the same length")
        }
    }

    if (!is.null(seq2)) {
        args1 <- paste(options, "--file-ref -x ", index, " -1 ", seq1, " -2 ", seq2,
            " -S ", output_file)
    } else {
        args1 <- paste(options, "--file-ref -x ", index, " -U ", seq1, " -S ", output_file)
    }

    # show(paste('nvBowtie ', args1))

    .callbinary2(bin1 = "nvBowtie", args1)
}



#' @name nvbio_tests
#' @title Print unit tests
#' @description Calling nvbio_tests() performs alignment tests
#' @author Franck RICHARD
#' @return The value is 0 when there is not any mistakes.
#' @export nvbio_tests
#' @examples
#' nvbio_tests()

nvbio_tests <- function() {
    if (R.Version()$arch == "i386") {
        return("nvbio_tests is not available for 32bit, please use 64bit R instead")
    }

    .callbinary2(bin1 = "nvbio-test", args1 = "-aln wfa")
}

#' @name nvBowtie_version
#' @title Print version information
#' @description Calling nvBowtie_version() displays the version information
#' of the RbowtieCuda package in use.
#' @author Franck RICHARD
#' @return An invisible \code{Integer} of call status.
#' The value is 0 when there is not any mistakes.
#' @export nvBowtie_version
#' @examples
#' nvBowtie_version()

nvBowtie_version <- function() {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit, please use 64bit R instead")
    }

    .callbinary2(bin1 = "nvBowtie", args1 = "--version")
}

#' @name nvBowtie_usage
#' @title Print available arguments that can be passed to nvBowtie()
#' @description Calling nvBowtie_usage() prints available arguments 
#' that can be passed to the nvBowtie() function.
#' @author Franck RICHARD
#' @return An invisible Integer of call status. The value is 0 when there is not any mistakes
#' Otherwise the value is non-zero.
#' @references Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read
#' alignment with Bowtie 2. Nature methods, 9(4), 357-359.
#' @export nvBowtie_usage
#' @importFrom methods show
#' @examples
#' nvBowtie_usage()



nvBowtie_usage <- function() {
    if (R.Version()$arch == "i386") {
        return("nvBowtie is not available for 32bit, please use 64bit R instead")
    }

    # invisible(.callbinary2( bin1 = 'nvBowtie', args1 = '-h'))

    show("options:")
    show("General:")
    show("  --verbosity         int                    [5]        verbosity level")
    show("  --upto       | -u   int                    [-1]       maximum number of reads to process")
    show("  --trim3      | -3   int                    [0]        trim the first N bases of 3'")
    show("  --trim5      | -5   int                    [0]        trim the first N bases of 5'")
    show("  --nofw                                     [false]    do not align the forward strand")
    show("  --norc                                     [false]    do not align the reverse-complemented strand")
    show("  --device            int                    [0]        select the given cuda device(s) (e.g. --device 0 --device 1 ...)")
    show("  --file-ref                                 [false]    load reference from file")
    show("  --server-ref                               [false]    load reference from server")
    show("  --phred33                                  [true]     qualities are ASCII characters equal to Phred quality + 33")
    show("  --phred64                                  [false]    qualities are ASCII characters equal to Phred quality + 64")
    show("  --solexa-quals                             [false]    qualities are in the Solexa format")
    show("  --rg-id             string                            add the RG-ID field of the SAM output header")
    show("  --rg                string,val                        add an RG-TAG field of the SAM output header")
    show("  --cache-writes      bool                   [false]    speed up writes on disk")
    show("Paired-End:")
    show("  --ff                                       [false]    paired mates are forward-forward")
    show("  --fr                                       [true]     paired mates are forward-reverse")
    show("  --rf                                       [false]    paired mates are reverse-forwardd")
    show("  --rr                                       [false]    paired mates are reverse-reverse")
    show("  --minins     |  -I  int                    [0]        minimum insert length")
    show("  --maxins     |  -X  int                    [500]      maximum insert length")
    show("  --overlap                                  [true]     allow overlapping mates")
    # show(' --no-discordant [false] do not allow discordant mates')
    show("  --no-mixed                                 [false]    only report paired alignments")
    show("  --ungapped-mates | -ug                                perform ungapped mate alignment")
    show("Seeding:")
    show("  --seed-len   | -L   int                    [22]       seed lengths")
    show("  --seed-freq  | -i   {G|L|S},x,y                       seed interval, as x + y*func(read-len) (G=log,L=linear,S=sqrt)")
    show("  --max-hits          int                    [100]      maximum amount of seed hits")
    show("  --max-reseed | -R   int                    [2]        number of reseeding rounds")
    show("  --top               bool                   [false]    explore top seed entirely")
    show("  --N                 bool                   [false]    allow substitution in seed")
    show("Extension:")
    show("  --mode              {best,best-exact,all}  [best]     alignment mode\n")
    show("  --all        | -a                          [false]    perform all-mapping (i.e. find and report all alignments)")
    show("  --local                                    [false]    perform local alignment")
    show("  --rand                                     [true]     randomized seed hit selection")
    show("  --no-rand                                  [false]    do not randomize seed hit selection")
    show("  --max-dist          int                    [15]       maximum edit distance")
    show("  --max-effort-init   int                    [15]       initial maximum number of consecutive extension failures")
    show("  --max-effort | -D   int                    [15]       maximum number of consecutive extension failures")
    show("  --min-ext           int                    [30]       minimum number of extensions per read")
    show("  --max-ext           int                    [400]      maximum number of extensions per read")
    show("  --very-fast                                           apply the very-fast presets")
    show("  --fast                                                apply the fast presets")
    show("  --sensitive                                           apply the sensitive presets")
    show("  --very-sensitive                                      apply the very-sensitive presets")
    show("  --very-fast-local                                     apply the very-fast presets")
    show("  --fast-local                                          apply the fast presets")
    show("  --sensitive-local                                     apply the sensitive presets")
    show("  --very-sensitive-local                                apply the very-sensitive presets")
    show("Scoring:")
    show("  --scoring           {sw|ed|wfa}            [ed]       Smith-Waterman / Edit-Distance / Wfa scoring")
    show("  --score-min         {G|L|S},x,y                       minimum score function, as x + y*func(read-len)")
    show("  --ma                int                               match bonus")
    show("  --mp                int,int                           mismatch min/max penalties")
    show("  --np                int                               N penalty")
    show("  --rdg               int,int                           read open/extension gap penalties")
    show("  --rfg               int,int                           reference open/extension gap penalties")
    show("Alternative:")
    show("  --wfa                                                 Activate wavefront algorithm")
    show("Reporting:")
    show("  --mapQ-filter | -Q  int                    [0]        minimum mapQ threshold")
    show("")
    show("")
    show("Default values are indicated in brackets [].")
    show("")
    show("The '--scoring-scheme filename' option allows to provide a custom Smith-Waterman scoring")
    show("scheme through a text file, where each line must contain a token value pair.")
    show(" The tokens and default values are reported below:")
    show("*  match               0        // local alignment: 2")
    show("*  mm-penalty-min      2")
    show("*  mm-penalty-max      6")
    show("*  N-penalty-min       1")
    show("*  N-penalty-max       1")
    show("*  score-min-const     -0.6     // local alignment: 0")
    show("*  score-min-coeff     -0.6     // local alignment: 10")
    show("*  score-min-type      linear   // local alignment: log")
    show("*  N-ceil-const        0")
    show("*  N-ceil-coeff        0.15")
    show("*  read-gap-const      5")
    show("*  read-gap-coeff      3")
    show("*  ref-gap-const       5")
    show("*  ref-gap-coeff       3")
    show("*  gap-free            5")
}



#' @name .callbinary2
#' @title Make system call for binaries
#' @description Function that makes a system call for the RbowtieCuda binaries.
#'  Note that it is not designed to be used outside the package.
#' @author Franck RICHARD
#' @param bin1 \code{Character}. The binary used for the system call.
#' @param args1 \code{Character}. The arguments to pass to the binary.
#' @param path \code{Character} Optional: If passed to function, returns
#'  the path.
#' @export .callbinary2
#' @return An invisible \code{Integer} of call status.
#' The value is 0 when there is not any mistakes.
#' @examples
#' .callbinary2(bin1 = "nvBowtie", args1 = "--version")

.callbinary2 <- function(bin1, args1, path = NULL) {

    bin1 <- paste(file.path(system.file(package = "RbowtieCuda"), bin1))

    result <- system2(bin1, args1)

    if (!is.null(path))
        return(path) else return(result)
}