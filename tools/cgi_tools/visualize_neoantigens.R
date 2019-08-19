library(ggplot2)
sink(stdout(), type = "message") # sink messages to stdout
require(grid)
library(extrafont)
sink(NULL, type="message") # close the sink
#font_import()

args <- commandArgs(trailing=TRUE)

neoantigen.file <- args[1]
pdf.file <- args[2]

neoantigens <- read.table(neoantigen.file, sep="\t", header=TRUE)
neo.filtered <- neoantigens[neoantigens $NM_FOLD_CHANGE > 1.5,]

p6 <- ggplot(neo.filtered, aes(x = neo.filtered$VAF , y = neo.filtered$AD, size = neo.filtered$NM_FOLD_CHANGE, fill = neo.filtered$FPKM)) + geom_point(shape = 21) + scale_size(range = c(1, 10)) + scale_fill_continuous(low = "blue", high = "red") + ggtitle("Neoantigen Ranking") + theme(plot.title = element_text(family="ArialMT", hjust = 0.5)) + labs(x = "Exome Variant Allele Frequency", y = "Allele Depth") + labs(size = "Binding Affinity Fold Change", fill = "fpkm") + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", legend.key.size = unit(1, "cm"), axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 14, family = "ArialMT", face = "bold"), text=element_text(family="ArialMT"), axis.text.x=element_text(colour="black", size = 9, family="ArialMT"), axis.text.y=element_text(family="ArialMT", colour="black", size = 9)) + geom_vline(xintercept=c(.40), linetype="dashed", color="red", show.legend=TRUE) + geom_hline(yintercept=c(60), linetype="dashed", color="darkgreen", show.legend=TRUE) + annotate("text", max(neo.filtered$VAF), 60, hjust = 1, vjust = -1, label = "Depth Cutoff") + annotate("text", 0.40, max(neo.filtered$AD), hjust = 1, vjust = -1, label = "Exome VAF Cutoff", angle=90) + geom_text(aes(label=ifelse(AD>=60 & VAF>=0.40,as.character(paste(SEQ_NAME,"(", format(round(FPKM, 2)), ",", format(round(NM_FOLD_CHANGE, 1)), ")")),'')),hjust=0,vjust=0, size=2)

#ggsave(pdf.file, plot=p6,  
#       width=5, height=3.5, units="in", device = pdf)

pdf(pdf.file)
p6+theme(text=element_text(family="ArialMT"))
dev.off()
