### want CNV data filtered into oncogenes and then take only amplifations
### want CNV data filtered into tumour suppressor genes and take only deletion/loss

setwd("~/Documents/R/Secrier_Data")

#fist step is to take our COSMIC list and extract oncogenes and tumour suppressor genes
cancer_genes_df <- read.csv("COSMICCancerGeneCensus_allMon_Oct_8_13_35_11 2018.csv")

##filtering for all cases that include either TSG or oncogene
can_gene_onco_all <- cancer_genes_df[which(grepl("oncogene",cancer_genes_df$Role.in.Cancer)),]
can_gene_TSG_all <- cancer_genes_df[which(grepl("TSG",cancer_genes_df$Role.in.Cancer)),]

save(can_gene_onco_all, file = "COSMIC_oncogenes.Rda")
save(can_gene_TSG_all, file = "COSMIC_TSG.Rda")


