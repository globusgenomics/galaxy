#!/usr/bin/env Rscript

library("optparse")
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--factors"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-o", "--outputF"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library("rjson")
parser<-newJSONParser()
parser$addData(opt$factors)
factorList <- parser$getObject()
#print(factorList)
suppressPackageStartupMessages({
  library('sleuth')
  library('wasabi')
})

# set the number of available cores to 4
options(mc.cores = 8L)

sf_dir <-dir(opt$dir)
sf_dir1<-paste(opt$dir,'/',sf_dir,sep='')
prepare_fish_for_sleuth(sf_dir1)
#print(sf_dir)

tmp<-summary(factorList)
condition <-NULL
for (i in range(1,nrow(tmp))){
  condition<-c(condition,rep(rownames(tmp)[i],tmp[i]))
}

#tmp<-names(factorList)
s2c<-data.frame(path=sf_dir1, sample=sf_dir, condition)
s2c$path <- as.character(s2c$path)

#print(s2c)

so <- sleuth_prep(s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
#models(so)

#for (i in 2:length(unique(condition))){
#  x<-paste('condition',unique(condition)[i],sep="")
#  so <- sleuth_wt(so, x)
#}

#so <- sleuth_wt(so, 'conditioncond2')
#results_table <- sleuth_results(so, 'conditioncond2'
#print(opt$outputF)
final_output <- results_table[order(results_table$qval),]
write.table(final_output,file=opt$outputF,row.names=F,sep = "\t",quote=FALSE)
