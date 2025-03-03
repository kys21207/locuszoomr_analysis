if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

install.packages("locuszoomr")

library(locuszoomr)
library(data.table)

input_file="L12_PSORIASIS_meta_out_filtered.tsv.gz"
# Read the gzipped input file
data <- fread(cmd = paste("gunzip -c", shQuote(input_file)), sep = "\t", header = TRUE)
        
MVP_EUR <- data[, .(chrom = `#CHR`, pos = POS, other_allele = REF, effect_allele = ALT, 
                    beta = MVP_EUR_beta, se = MVP_EUR_sebeta, p = MVP_EUR_pval, 
                    r2 = MVP_EUR_r2, rsid = rsid)]
MVP_EUR <- na.omit(MVP_EUR)           

MVP_AFR <- data[, .(chrom = `#CHR`, pos = POS, other_allele = REF, effect_allele = ALT, 
                    beta = MVP_AFR_beta, se = MVP_AFR_sebeta, p = MVP_AFR_pval, 
                    r2 = MVP_AFR_r2, rsid = rsid)]
MVP_AFR <- na.omit(MVP_AFR)           

MVP_AMR <- data[, .(chrom = `#CHR`, pos = POS, other_allele = REF, effect_allele = ALT, 
                    beta = MVP_HIS_beta, se = MVP_HIS_sebeta, p = MVP_HIS_pval, 
                    r2 = MVP_HIS_r2, rsid = rsid)]
MVP_AMR <- na.omit(MVP_AMR)           

library(EnsDb.Hsapiens.v86)
loc1 <- locus(data = MVP_EUR, gene = 'TYK2', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86")
loc1 <- link_LD(loc1, pop = "CEU", genome_build="grch38", token = "abc8c8e6decc")
summary(loc1)
loc2 <- locus(data = MVP_AFR, gene = 'TYK2', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86")
loc2 <- link_LD(loc2, pop = "LWK", genome_build="grch38" ,token = "abc8c8e6decc")
summary(loc2)

loc3 <- locus(data = MVP_AMR, gene = 'TYK2', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86")
loc3 <- link_LD(loc3, pop = "MXL", genome_build="grch38" ,token = "abc8c8e6decc")
summary(loc3)


pdf("TYK2_plot_MVP.pdf", width = 6, height = 8)
# set up layered plot with 2 plots & a gene track; store old par() settings
oldpar <- set_layers(3)
scatter_plot(loc1, xticks = FALSE, labels="index", main="MVP EUR", ylim = c(0, 35))
scatter_plot(loc2, xticks = FALSE, labels="index", main="MVP AFR", ylim = c(0, 35))
scatter_plot(loc3, xticks = FALSE, labels="index", main="MVP AMR", ylim = c(0, 35))
genetracks(loc1)
par(oldpar)  # revert par() settings
dev.off()
