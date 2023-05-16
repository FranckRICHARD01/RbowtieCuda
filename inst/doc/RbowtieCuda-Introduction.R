## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("RbowtieCuda")

## ----loading------------------------------------------------------------------
library(RbowtieCuda)

## ----idad---------------------------------------------------------------------
td <- tempdir()
fa_file <- system.file(package="RbowtieCuda", "extdata", "bt2", "refs", "lambda_virus.fa")
nvBWT(myinput=fa_file, output=file.path(td, "index"), options="")

## ----bt2bd1-------------------------------------------------------------------

read_1 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_1.fastq")
read_2 <- system.file(package="RbowtieCuda", "extdata", "bt2", "reads", "reads_2.fastq")
nvBowtie(file.path(td, "index"), file.path(td, "my_result.bam"), options="", seq1=read_1, seq2=read_2)

## ----bt2usage---------------------------------------------------------------
nvBowtie_usage()

## ----bt2version---------------------------------------------------------------
nvBowtie_version()

## ----sessioninfo--------------------------------------------------------------
sessionInfo()

