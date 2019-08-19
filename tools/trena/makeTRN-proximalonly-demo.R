#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(TReNA)

# 0. set options
#genome.db.uri = "postgres://whovian/hg38"
#project.db.uri = "postgres://whovian/lymphoblast"
#out = "/proj/price1/sament/lymphoblast_trn/hg38.lymphoblast"
output_dir = args[1]
#print(output_dir)
promoter_counts = local(get(load(args[2])))
expr = local(get(load(args[3])))
method = args[4]

# 1, get counts of binding sites for each TF proximal to each gene (by default +/- 10kb from the TSS)

# 2. build TRN model by integrating TFBS counts with expression data

trn = makeTrnFromPromoterCountsAndExpression( counts = promoter_counts , expr = expr ,  method = method)

save( trn , file = paste( output_dir, "/trn.RData" , sep = "" ))


