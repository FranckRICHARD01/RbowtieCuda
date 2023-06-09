% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bowtieCuda.R
\name{nvBowtie}
\alias{nvBowtie}
\title{nvBowtie align function}
\usage{
nvBowtie(index, output_file, options, seq1, seq2 = NULL)
}
\arguments{
\item{index}{\code{Character} vector. Index file created by nvBWT.}

\item{output_file}{\code{Character} vector. The alignment output file (.sam o .bam)}

\item{options}{\code{Character} vector. Indicate additional options here to modify the processing. The complete
list of available options can be accessed by typing \code{nvBowtie_usage()}.}

\item{seq1}{\code{Character} vector. For single-end sequencing, it contains sequence file paths.
For paired-end sequencing, it can be file path with #1 mates paired with file paths in seq2.}

\item{seq2}{\code{Character} vector. It contains file paths with #2 mates paired with file path in seq1.
For single-end sequencoing files, it must be \code{NULL}.'}
}
\value{
An invisible \code{Integer} of call
status. The value is 0 when there is not any mistakes
Otherwise the value is non-zero.
}
\description{
nvBowtie is a GPU-accelerated re-engineering of Bowtie2, a very widely used short-read aligner.
}
\examples{
td <- tempdir()

## Building index
fa_file <- system.file(package="RbowtieCuda", "extdata", "bt2", "refs", "lambda_virus.fa")
nvBWT(myinput=fa_file, output=file.path(td, "index"), options="")

## Alignments
read_1 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_1.fastq")
read_2 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_2.fastq")

## Sam file created with paired-end files (paired mates are here reverse-forward)
nvBowtie(file.path(td, "index"), file.path(td, "my_result.sam"), options="--rf", seq1=read_1, seq2=read_2)

## Bam file created with single-end file
nvBowtie(file.path(td, "index"), file.path(td, "my_result.bam"), options="", seq1=read_1)
}
\references{
Langmead, B., & Salzberg, S. L. (2012).
Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
}
