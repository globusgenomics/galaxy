rm(list=ls())
options(warn=-1)
args=commandArgs()[-c(1:4)]
infile=args[1]
outfile=sub("Temp.2R$", "var", infile)

x=read.table(infile, sep="\t", header=FALSE, row.names=NULL, colClasses=c('factor', rep('integer', 2), rep('factor', 3), rep('integer', 5), 'factor'),
	col.names=c('chr', 'start', 'end', 'ref', 'alt', 'zyg', 'QC', 'NRF', 'NRR', 'NAF', 'NAR', 'VCF.FILTER'))

x=cbind(x, matrix(
	as.numeric(format(sapply(1:nrow(x), function(i)
	with(x[i, ], c(with(binom.test(c(sum(NRF, NAF), sum(NRR, NAR))), p.value),
	with(binom.test(c(NRF, NRR)), p.value),
	with(binom.test(c(NAF, NAR)), p.value)))), digits=3, scientific=TRUE)),
	ncol=3, byrow=TRUE, dimnames=list(NULL, c('p.strand', 'p.ref', 'p.alt')))
)

x=cbind(subset((x=cbind(x, do.call('rbind', lapply(1:nrow(x), function(i) with(x[i, ], {
	AD=sum(NAF, NAR)
	DP=sum(NRF, NRR, AD)
	AF=signif(AD/DP, digits=3)
	VAR.FILTER=c(zyg!='unk' & QC>=20 & DP>=10 & AD>=5 & (p.strand>.05 | min(sum(NRF, NAF), sum(NRR, NAR))>=10) & (p.ref>.05 | min(NRF, NRR)>=10) & (p.alt>.05 | min(NAF, NAR)>=10))
	cbind(DP,AD,AF, VAR.FILTER)
}))))), select=-VAR.FILTER), VAR.FILTER=with(x, factor(VAR.FILTER, levels=c(0, 1), labels=c('SKIP', 'PASS'))))

write.table(as.matrix(x), outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
q(runLast=FALSE)