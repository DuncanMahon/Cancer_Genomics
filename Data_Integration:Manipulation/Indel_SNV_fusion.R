###Fusing INDEL and SNV data for later analysis

#load in indels data
setwd("~/Documents/R/Secrier_Data")
load("indels_129cohort.RData")

#adding shortened illumina ID to indel data
indels$sample <- sapply(indels$V2,
                      function(x) strsplit(x,"_vs_")[[1]][1])
colnames(indels)[11] <- "type"

#reading in cancer genes for filtering
cancer_genes_df <- read.csv("COSMICCancerGeneCensus_allMon_Oct_8_13_35_11 2018.csv")
cancer_genes <- cancer_genes_df[,1]

#filtering for cancer genes in indel data
indels_cg <- indels[which(indels$gene %in% cancer_genes),]

#load SNV data
load("snvs_129cohort.RData")
#we just want missense SNV's
snvs.filtered <- snvs[which(grepl("missense",snvs$Consequence)),]
snvs.missense <- unique(snvs[,c("Sample_ID","sample","gene")])
snvs.missense <- snvs.missense[which(snvs.missense$gene != "-"),]
snvs.selected <- snvs.missense[which(snvs.missense$gene %in% cancer_genes),]

#filtered indel data by cancer genes
#first cut down our indel data
indels_cutdown <- data.frame("sample" = indels$sample,"gene" = indels$gene,"type" = indels$type)
#strategy: grepl each type seperately then bind the 3 new dataframes
indels_fv <- indels_cutdown[which(grepl("frameshift_variant",indels_cutdown$type)),]
indels_id <- indels_cutdown[which(grepl("inframe_deletion",indels_cutdown$type)),]
indels_ii <- indels_cutdown[which(grepl("inframe_insertion",indels_cutdown$type)),]
indels_filtered <- rbind(indels_fv,indels_id,indels_ii)

#merge the SNV and the indel data
#add type column to snv dataframe and remove Sample_ID
snvs.missense$type <- "SNV"
snvs.missense$Sample_ID <- NULL

#bind the two together
SNV_Indel_data <- rbind(snvs.missense,indels_filtered)

## refining by removing type and doing uniques instances for the whole df
SNV_Indel_data <- within(SNV_Indel_data, rm(type))
SNV_Indel_data <- unique(SNV_Indel_data)
#save SNV_Indel_data
save(SNV_Indel_data, file="SNV_Indel_Combined_COSMIC.Rda")



