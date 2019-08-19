#!/usr/bin/env Rscript

suppressMessages(library(optparse))

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("", "--fwd"), dest="fnFs",
                help="The path(s) to the input fastq file(s)."),
    make_option(c("", "--rev"), dest="fnRs",
                help="The path(s) to the input reverse fastq file(s) from paired-end sequence data corresponding to those provided to the fwd argument. If NULL, the fwd files are processed as single-reads."),
    make_option(c("", "--filt"), dest="filtFs",
                help="The path(s) to the output filtered file(s) corresponding to the fwd input files. If containing directory does not exist, it will be created."),
    make_option(c("", "--filt.rev"), dest="filtRs",
                help="The path(s) to the output fastq file(s) corresponding to the rev input. Can also provide a directory, which if not existing will be created (how to differentiate between dir/file in len1 case?)."),
    make_option(c("", "--truncQ"), dest="truncQ", default=2,
                help="Truncate reads at the first instance of a quality score less than or equal to truncQ."),
    make_option(c("", "--truncLen"), dest="truncLen", default=0,
                help="Truncate reads after truncLen bases. Reads shorter than this are discarded."),
    make_option(c("", "--trimLeft"), dest="trimLeft", default=0,
                help="The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft."),
    make_option(c("", "--maxLen"), dest="maxLen", default=Inf,
                help="Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation."),
    make_option(c("", "--minLen"), dest="minLen", default=20,
                help="Remove reads with length less than minLen. minLen is enforced after trimming and truncation."),
    make_option(c("", "--maxN"), dest="maxN", default=0,
                help="After truncation, sequences with more than maxN Ns will be discarded. Note that dada does not allow Ns."),
    make_option(c("", "--minQ"), dest="minQ", default=0,
                help="After truncation, reads contain a quality score less than minQ will be discarded."),
    make_option(c("", "--maxEE"), dest="maxEE", default=Inf,
                help="After truncation, reads with higher than maxEE \"expected errors\" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))"),
    make_option(c("", "--primer.fwd"), dest="primer.fwd", default=NULL,
                help="Paired-read filtering only. A character string defining the forward primer. Only allows unambiguous nucleotides. The primer will be compared to the first len(primer.fwd) nucleotides at the start of the read. If there is not an exact match, the read is filtered out. For paired reads, the reverse read is also interrogated, and if the primer is detected on the reverse read, the forward/reverse reads are swapped."),
    make_option(c("", "--matchIDs"), dest="matchIDs", action="store_true", default=FALSE,
                help="Paired-read filtering only. Whether to enforce matching between the id-line sequence identifiers of the forward and reverse fastq files. If TRUE, only paired reads that share id fields (see below) are output. If FALSE, no read ID checking is done. Note: matchIDs=FALSE essentially assumes matching order between forward and reverse reads. If that matched order is not present future processing steps may break (in particular mergePairs)."),
    make_option(c("", "--id.sep"), dest="id.sep", default="\\s",
                help="Paired-read filtering only. The separator between fields in the id-line of the input fastq files. Passed to the strsplit."),
    make_option(c("", "--id.field"), dest="id.field", default=NULL,
                help="Paired-read filtering only. The field of the id-line containing the sequence identifier. If NULL (the default) and matchIDs is TRUE, the function attempts to automatically detect the sequence identifier field under the assumption of Illumina formatted output."),
    make_option(c("", "--multithreaded"), dest="multithreaded", action="store_true", default=FALSE,
                help="If TRUE, input files are filtered in parallel via mclapply. If an integer is provided, it is passed to the mc.cores argument of mclapply. Note that the parallelization here is by forking, and each process is loading another fastq file into memory. This option is ignored in Windows, as Windows does not support forking, with mc.cores set to 1. If memory is an issue, execute in a clean environment and reduce the chunk size n and/or the number of threads."),
    make_option(c("", "--n"), dest="n", default=1e5,
                help="The number of records (reads) to read in and filter at any one time. This controls the peak memory requirement so that very large fastq files are supported. See FastqStreamer for details."),
    make_option(c("", "--verbose"), dest="verbose", action="store_true", default=FALSE,
                help="Whether to output status messages.")
)

opt <- parse_args(OptionParser(option_list=option_list), args)

opt$fnFs <- strsplit(as.character(opt$fnFs), ',')[[1]]
opt$fnRs <- strsplit(as.character(opt$fnRs), ',')[[1]]
if(!is.null(opt$filtFs)){
    opt$filtFs <- strsplit(as.character(opt$filtFs), ',')[[1]]
}
if(!is.null(opt$filtRs)){
    opt$filtRs <- strsplit(as.character(opt$filtRs), ',')[[1]]
}
print(opt)

suppressMessages(library(dada2))

filterAndTrim(
    opt$fnFs,
    opt$filtFs,
    opt$fnRs,
    opt$filtRs,
    truncLen=opt$truncLen,
    maxN=opt$maxN,
    maxEE=opt$maxEE,
    truncQ=opt$truncQ,
    rm.phix=TRUE,
    compress=TRUE,
    multithread=opt$multithreaded
)
