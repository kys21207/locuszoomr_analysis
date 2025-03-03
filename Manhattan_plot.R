install.packages("qqman")
library(qqman)

png("manhattan_plot.png", width = 1000, height = 600)
 manhattan(MVP_AFR, chr="chrom", bp="pos", snp="rsid", p="p" )
dev.off()
