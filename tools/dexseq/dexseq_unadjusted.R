library(DEXSeq)

cargs <- commandArgs()
cargs <- cargs[(which(cargs == "--args")+1):length(cargs)]

annotationfile <- cargs[1]
covariates_input <- cargs[2] # Condition,RIN,numeric.batch
covariates <- noquote(strsplit(covariates_input, ","))[[1]]

formula.covs_input <- cargs[3] # RIN,numeric.batch
formula.covs <- noquote(strsplit(formula.covs_input, ","))[[1]]

filenames_and_samplenames <- cargs[6:length(cargs)-2]
filenames <- filenames_and_samplenames[seq(1, length(filenames_and_samplenames), 3)]
sample_names <- filenames_and_samplenames[seq(2, length(filenames_and_samplenames), 3)]
conditions_input <- filenames_and_samplenames[seq(3, length(filenames_and_samplenames), 3)]
conditions <- noquote(strsplit(conditions_input, ","))

run.dexseq.unadjusted.all.genes <- function() {

  # covariates file for the PD and control samples that passed QC
  samples.factors <- data.frame(row.names=sample_names, matrix(unlist(conditions), nrow=length(conditions), byrow=T))
  colnames(samples.factors) <- covariates

  #setwd("/unprotected/projects/mlpd/Andra/DEXSeq/input/htseq-count_exons")

  # create the ExonCountSet object
  print("reading ExonCountSet object")
  ecs <- read.HTSeqCounts(countfiles = filenames, design = samples.factors, flattenedfile = annotationfile)
  print("read ExonCountSet object")  

  # estimate size factors
  ecs <- estimateSizeFactors(ecs)
  print("estimated size factors")

  # estimate dispersions
  print("estimating dispersions")
  ecs <- estimateDispersions(ecs, nCores=8)
  print("estimated dispersions")

  # fit dispersions
  ecs <- fitDispersionFunction(ecs)
  print("fitted dispersions")

  # plot dispersion estimates
  pdf("all.genes.DispEsts.pdf")
  meanvalues <- rowMeans(counts(ecs))
  plot(meanvalues, fData(ecs)$dispBeforeSharing, log = "xy", main = "mean vs CR dispersion", xlab = "mean of normalized counts", ylab = "dispersion")
  x <- 0.01:max(meanvalues)
  y <- ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2]/x 
  lines(x, y, col = "red")
  dev.off()

  # test for differential expression
  print("testing for DEU")
  ecs <- testForDEU(ecs, nCores=8)
  print("tested for DEU")
  ecs <- estimatelog2FoldChanges(ecs)
  results <- DEUresultTable(ecs)
  write.table(results, "DEU_results_all_genes.txt", sep="\t", quote=F)

  pdf("all.genes.DEU.results.pdf")
  plot(results$meanBase, results[, "log2fold(PD/Control)"], log = "x", col = ifelse(results$padjust < 0.05, "red", "black"), ylim = c(-4, 4), main = "PD/C all genes' exons MA plot", ylab = "log2fold(PD/C)", xlab = "mean of normalized counts")
  dev.off()
  
  print("working on HTML")
  DEXSeqHTML(ecs, path="./",  FDR=0.1, color=c("#FF000080", "#0000FF80"))
  print("HTML done")
  #save.image("all_genes_unadjusted/DEXSeq.all.genes.exons.RData")

  return(ecs)
}


resulting.deu.object <- run.dexseq.unadjusted.all.genes()

print("finished run - starting to save the RData file")
save.image("DEXSeq.all.genes.exons.RData")
print("RData file saved")
