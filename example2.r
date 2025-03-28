if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

install.packages("locuszoomr")

library(locuszoomr)
library(data.table)
input_file1="/mnt/project/approved_results/pd_w_age_bsize_400_firth_parkinsons_all_comers_v3.regenie.tsv.gz"
data1 <- fread(cmd = paste("gunzip -c", shQuote(input_file1)), sep = " ", header = TRUE)

input_file2="/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Ast.tsv.gz"
data2 <- fread(cmd = paste("gunzip -c", shQuote(input_file2)), sep = "\t", header = TRUE)

input_file3="/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Mic.tsv.gz"
data3 <- fread(cmd = paste("gunzip -c", shQuote(input_file3)), sep = "\t", header = TRUE)

input_file4="/mnt/project/publically_available_supporting_files/pops_data/demo_file_rsID4regenie.gwaslab.gz"
ref <- fread(cmd = paste("gunzip -c", shQuote(input_file4)), sep = "\t", header = TRUE)

merged_data1 <- merge(data1, ref, by.x ="ID",by.y="SNPID")
merged_data1 <- merged_data1[!is.na(rsID)]
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO
PD <- merged_data1[, .(chrom = CHROM, pos = GENPOS, other_allele = ALLELE0, effect_allele = ALLELE1, 
                    beta = BETA, se = SE, p = 10^(-LOG10P), rsid = rsID)]

#celltype        gene_symbol     gene_id snps    chr38   pos38   REF     ALT     ALT_AF  beta    se      pvalue  significant_by_2step_FDR
AST <- data2[gene_id == "ENSG00000157601", .(chrom =  gsub("chr", "", chr38), pos = pos38, other_allele = REF, effect_allele = ALT, 
                    beta=beta,se=se,p = pvalue, rsid = snps)]
MIC <- data3[gene_id == "ENSG00000157601", .(chrom =  gsub("chr", "", chr38), pos = pos38, other_allele = REF, effect_allele = ALT, 
                    beta=beta,se=se,p = pvalue, rsid = snps)]
                    library(EnsDb.Hsapiens.v86)
loc1 <- locus(data = PD, gene = 'MX1', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86")
loc1 <- link_LD(loc1, pop = "CEU", genome_build="grch38", token = "abc8c8e6decc")
summary(loc1)
loc2 <- locus(data = AST, gene = 'MX1', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86")
loc2 <- link_LD(loc2, pop = "CEU", genome_build="grch38" ,token = "abc8c8e6decc")
summary(loc2)


pdf("MX1_plot_PD_vs_Ast.pdf", width = 4, height = 6)
# set up layered plot with 2 plots & a gene track; store old par() settings
oldpar <- set_layers(2)
scatter_plot(loc1, xticks = FALSE, labels="index", main="Parkinson Disease", ylim = c(0, 10), beta = "beta")
scatter_plot(loc2, xticks = FALSE, labels="rs17000900", main="Astrocytes.", ylim = c(0, 10), beta = "beta")
genetracks(loc1)
par(oldpar)  # revert par() settings
dev.off()
