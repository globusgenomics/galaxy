#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(TReNA)
output_dir = args[1]
genome.db.uri = args[2]
project.db.uri = args[3]
ncores=8
# 1, get counts of binding sites for each TF proximal to each gene (by default +/- 10kb from the TSS)


promoter_counts = getTfbsCountsInPromoters( 
   genome.db.uri=genome.db.uri , project.db.uri=project.db.uri , 
   size.upstream = 10000 , size.downstream = 10000, cores = ncores  # define the size of the window around the TSS
) # specify the number of cores for parallelization
save(promoter_counts , file = paste(output_dir,"/promoter_tfbs_counts.RData",sep=""))
