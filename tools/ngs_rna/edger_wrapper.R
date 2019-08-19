#!/usr/bin/Rscript

suppressPackageStartupMessages(library("edgeR")) # Ignore "Loading required package" non error statement that is printed to STDERR.   (This makes Galaxy happy.)
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("Cairo"))

#  1.  Read in list of files and conditions
args <- commandArgs(trailingOnly = TRUE)       #Reads the command line arguments in order to find the file listing the count tables and conditions

list <- args[1]                    #Assigns the list to an object
mds_plot_fn <- args[2]
tagwise_disp_plot_fn <- args[3]
ma_plot_fn <- args[4]
diff_expr_fn <- args[5]
sig_diff_expr_fn <- args[6]


targets <- read.delim(file = list, stringsAsFactors = FALSE)   #Reads in the list and assigns it to the targets object

#  2.  Create the DGEList object for differential testing
d <- readDGE(targets)     #Calculates the size of the count libraries and produces the DGEList object
cpm.d <- cpm(d)                                                                 #  Compute counts per million
samples = d$samples$description                                                 #  Get sample names to use as labels below
colnames(d) <- c(samples)                                                       #Assigns the samplenames to the column headers for further analysis

#  3.  Filter out genes that do not reach a minimum CPM threshold.
MIN_REP_NUMBER  = min(table(d$samples$group))                           #  Determine mimimum number of replicates for a condition
CPM_THRESHOLD   = "5"                           #  Set minimum CPM threshold
EXPRESSED = rowSums(cpm.d >= CPM_THRESHOLD) >= MIN_REP_NUMBER   #  Find genes that reach CPM_TRESHOLD in at least MIN_REP_NUMBER samples
d <- d[EXPRESSED, ]                                             #  Keep genes that meet the threshold
dim(d)

#  4.  Calculate effective library size
d <- calcNormFactors(d)         #  There are multiple ways to correct for differences in RNA composition

#  5.  Make an MDS plot
options(bitmapType='cairo')
#png(file=mds_plot_fn, height=600, width=600)
Cairo(file=mds_plot_fn, height=600, width=600, type="png")
MDS = plotMDS(d, main="MDS Plot", labels=samples)              #  Make an MDS plot.  Save the object and then make a formatted version
plot(MDS$x,MDS$y,pch=18,cex=2, main="MDS Plot",xlab='Dimension 1',ylab='Dimension 2',xlim=c(min(MDS$x)-0.25,max(MDS$x)+0.25),ylim=c(min(MDS$y)-1,max(MDS$y)+1));text(MDS$x,MDS$y+0.05,labels=samples,pos=3,cex=1);dev.off()

#  6.  Calculate and plot the tagwise dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)                                     #  Various parameters are available to control how strongly to weight the gene's dispersion vs. the dispersion trend
#png(file=tagwise_disp_plot_fn, height=600, width=600);plotBCV(d);dev.off()
Cairo(file=tagwise_disp_plot_fn, height=600, width=600,type="png");plotBCV(d);dev.off()

#  7.  Test for differential expression
DE <- exactTest(d)
n_GENES_TESTED = dim(DE$table)[1]
DIFF = topTags(DE,n=n_GENES_TESTED)$table               #  Get the differential expression results for all tested genes

#  8.  Make an MA plot, showing genes that meet the FDR threshold in red
FDR_THRESHOLD = 0.05
#png(file=ma_plot_fn, height=600, width=600);par(cex.lab=1.5,cex.axis=1.5,mar=c(5,5,5,5));plot(DIFF$logCPM,DIFF$logFC,xlab='log2(CPM)',ylab="logFC",pch=18,cex=1,col=ifelse(DIFF$FDR < FDR_THRESHOLD,'red','black'));abline(h = c(-2, 2), col = "dodgerblue");dev.off()
Cairo(file=ma_plot_fn, height=600, width=600, type="png");par(cex.lab=1.5,cex.axis=1.5,mar=c(5,5,5,5));plot(DIFF$logCPM,DIFF$logFC,xlab='log2(CPM)',ylab="logFC",pch=18,cex=1,col=ifelse(DIFF$FDR < FDR_THRESHOLD,'red','black'));abline(h = c(-2, 2), col = "dodgerblue");dev.off()

#  9.  Filter for significant results
DIFF_SIG <-  DIFF[DIFF$FDR < 0.05, ]
dim(DIFF_SIG)

#  10.  Combine data and write to file
#  -  d$pseudo.counts:          contains pseudo counts for each sample, in counts per million after correcting for library size
#  -  DIFF:             contains differential expression results for all tested genes
NORM_CPM = d$pseudo.counts
OUT_DATA = merge(NORM_CPM,DIFF,by.x="row.names",by.y="row.names")
names(OUT_DATA)[1] = "Gene_ID"
OUT_DATA = format(OUT_DATA, digits=5)
write.table(OUT_DATA, file=diff_expr_fn,row.names=TRUE, sep = "\t", quote = FALSE)

OUT_SIG = merge(NORM_CPM, DIFF_SIG, by.x="row.names",by.y="row.names")
names(OUT_SIG)[1] = "Gene_ID"
OUT_SIG = format(OUT_SIG, digits=5)
write.table(OUT_SIG, file=sig_diff_expr_fn, row.name=TRUE, sep = "\t", quote = FALSE)
