args <- commandArgs(trailingOnly = TRUE)
input_bed <- args[1]
input_vcf <- args[2]
output_vcf <- args[3]

x <- read.csv2( input_bed, sep='\t', header=FALSE)
y <- read.csv2( input_vcf, sep='\t', header=FALSE)

colnames(x)[c(1,3)] <- c('chr', 'loc')
colnames(y)[c(1,2)] <- c('chr', 'loc')

extract <- y[match(x$V3, y$V2),]

extract <- merge(x,y,by=c('chr','loc'))
extract <- extract[,c(1,2,8:15)]

write.table(extract, file=output_vcf, sep='\t', quote=FALSE, col.names=F, row.names=F)
