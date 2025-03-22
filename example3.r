if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

install.packages("locuszoomr")

library(locuszoomr)
library(data.table)
library(EnsDb.Hsapiens.v86)

source("scatter_plot_with_cs.R")

input_file1="/mnt/project/publically_available_supporting_files/gwas_public_results/meta_mvp_ukb_psoriasis.regenie.tsv.gz"
data1 <- fread(cmd = paste("gunzip -c", shQuote(input_file1)), sep = "\t", header = TRUE)

#input_file2="/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Ast.tsv.gz"
#data2 <- fread(cmd = paste("gunzip -c", shQuote(input_file2)), sep = "\t", header = TRUE)

#input_file3="/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Mic.tsv.gz"
#data3 <- fread(cmd = paste("gunzip -c", shQuote(input_file3)), sep = "\t", header = TRUE)

input_file4="/mnt/project/publically_available_supporting_files/pops_data/demo_file_rsID4regenie.gwaslab.gz"
ref <- fread(cmd = paste("gunzip -c", shQuote(input_file4)), sep = "\t", header = TRUE)

input_file5="susie_results/meta_mvp_ukb_psoriasis_credible_sets_block_35.txt"
cs <- read.table(input_file5,header=T,sep="\t")

merged_data1 <- merge(data1, ref, by.x ="ID",by.y="SNPID")
merged_data1 <- merged_data1[!is.na(rsID)]
merged_data2 <- merge(merged_data1, cs, by.x="ID", by.y="SNP",all.x=T)

#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE CHISQ LOG10P INFO
dta <- merged_data2[, .(chrom = CHROM, pos = GENPOS, other_allele = ALLELE0, effect_allele = ALLELE1, 
                    beta = BETA, se = SE, p = 10^(-LOG10P), rsid = rsID,cs_num = CS_Number)]
loc1 <- locus(data = dta, gene = 'IL23R', flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v86") 
pf <- quote({
  v <- loc1$TX[loc1$TX$gene_name == "IL23R", c("start", "end")]
  abline(v = v, col = "green")
})
oldpar <- set_layers(1)
scatter_plot_cs(loc1,panel.first = pf)
genetracks(loc1, highlight = "IL23R")
par(oldpar)
