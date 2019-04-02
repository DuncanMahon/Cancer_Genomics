### Exploring enrichment of SNV's on each CR

#remember to setwd
setwd("~/Documents/R/Secrier_Data")

load("snvs_129cohort.RData")
colnames(snvs)[1:11] <- c("OCCAMS_ID","Sample_ID","Mutation","Chromosome","Position",
                          "Ref","Alt","GeneID","TranscriptID","Type",
                          "Consequence")
save(snvs, file="snvs_129cohort.RData")

snvs.filtered <- snvs[which(grepl("missense",snvs$Consequence)),]
snvs.missense <- unique(snvs[,c("Sample_ID","sample","gene")])
snvs.missense <- snvs.missense[which(snvs.missense$gene != "-"),]

#reading clinical data
load("datOAC.clinicalPlusBFB.RData")

#reading  (and extracting) gene list
cancer_genes_df <- read.csv("COSMICCancerGeneCensus_allMon_Oct_8_13_35_11 2018.csv")
cancer_genes <- cancer_genes_df[,1]
snvs.selected <- snvs.missense[which(snvs.missense$gene %in% cancer_genes),]

#count frequency in these genes
table(snvs.selected$gene)
#length of this table gives 697 genes out of cancer gene list of 719

### Merge two tables:
dat.clinical.plusSNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                              "HasBFB","BFBcateg")],
                              snvs.selected[,c("sample","gene")],
                              by.x="Illumina_Barcode", by.y="sample",
                              all.x=TRUE, all.y=FALSE) 

### Count the frequency of mutations across BFB types:
tb <- table(dat.clinical.plusSNV[,c("BFBtype","gene")])

#all BFB types
BFB_types <- sort(unique(dat.clinical.plusSNV$BFBtype))
Gene_list <- variable.names(tb)

#function for pvalues over row from fisher test
multiple_fisher_thru_row <- function(t,n){
  p_vals <- c()
  for (i in 1:ncol(t)){
    x <- t[n,i]
    y <- rowSums(t)[n] - x
    z <- colSums(t)[i] - x
    xyz <- x + y + z
    br <- sum(t) - xyz
    m <- rbind(c(x,y),c(z,br))
    f <- fisher.test(m)
    f_pval <- f$p.value
    p_vals <- append(p_vals,f_pval)
  }
  return(p_vals)
}

BFB_like_pvals <- multiple_fisher_thru_row(t = tb, n = 1)
BFB_like_pvals_adjusted <- p.adjust(BFB_like_pvals,method = "BH")
df <- data.frame(Gene_list,BFB_like_pvals,BFB_like_pvals_adjusted)

#creating list of dataframes for all BFB_types
BFBtypes_df_list <- list()
for (i in 1:length(BFB_types)){
  pval <- multiple_fisher_thru_row(tb,i)
  pval_adjusted <- p.adjust(pval, method = "BH")
  BFBtypes_df_list[[i]] <- data.frame(Gene_list,pval,pval_adjusted)
}
BFBtypes_df_list
names(BFBtypes_df_list) <- BFB_types

#saving my tables
BFBlike_data <- BFBtypes_df_list[1]
save(BFBlike_data, file="BFBlike_analysis.Rda")
DoubleMinute_data <- BFBtypes_df_list[2]
save(DoubleMinute_data, file="DoubleMinute_analysis.Rda")
FocalAmplifications_data <- BFBtypes_df_list[3]
save(FocalAmplifications_data, file="FocalAmplifications_analysis.Rda")
No_BFB_data <- BFBtypes_df_list[4]
save(No_BFB_data, file="NO_BFB_analysis.Rda")
Subtelomeric_data <- BFBtypes_df_list[5]
save(Subtelomeric_data, file = "Subtelomeric_analysis.Rda")
TandemDuplication_data <- BFBtypes_df_list[6]
save(TandemDuplication_data, file = "TandemDuplication_analysis.Rda")

#analysis to find any significant relationships
#BFBlike_data$pval_adjusted

#loading in results
setwd("~/Documents/R/Secrier_Data")
load("BFBlike_analysis.Rda")
load("DoubleMinute_analysis.Rda")
load("FocalAmplifications_analysis.Rda")
load("NO_BFB_analysis.Rda")
load("Subtelomeric_analysis.Rda")
load("TandemDuplication_analysis.Rda")

#Extracting results
for (i in 1:length(BFBlike_data$`BFB-like`$pval_adjusted)){
  if (BFBlike_data$`BFB-like`$pval_adjusted[i] <= 0.05){
    print(BFBlike_data$`BFB-like`$Gene_list[i])
  }
}
#none
for (i in 1:length(DoubleMinute_data$DoubleMinute$pval_adjusted)){
  if (DoubleMinute_data$DoubleMinute$pval_adjusted[i] <= 0.05){
    print(DoubleMinute_data$DoubleMinute$Gene_list[i])
  }
}
#none
for (i in 1:length(FocalAmplifications_data$FocalAmplification$pval_adjusted)){
  if (FocalAmplifications_data$FocalAmplification$pval_adjusted[i] <= 0.05){
    print(FocalAmplifications_data$FocalAmplification$Gene_list[i])
  }
}
#none
for (i in 1:length(No_BFB_data$None$pval_adjusted)){
  if (No_BFB_data$None$pval_adjusted[i] <= 0.05){
    print(No_BFB_data$None$Gene_list[i])
  }
}
#none
for (i in 1:length(Subtelomeric_data$Subtelomeric$pval_adjusted)){
  if (Subtelomeric_data$Subtelomeric$pval_adjusted[i] <= 0.05){
    print(Subtelomeric_data$Subtelomeric$Gene_list[i])
  }
}
#none
for (i in 1:length(TandemDuplication_data$TandemDuplication$pval_adjusted)){
  if (TandemDuplication_data$TandemDuplication$pval_adjusted[i] <= 0.05){
    print(TandemDuplication_data$TandemDuplication$Gene_list[i])
  }
}
#none
