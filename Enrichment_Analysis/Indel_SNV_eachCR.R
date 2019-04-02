#Indel + SNV analysis (each BFB)

#load indel/snv data + clinical data
setwd("~/Documents/R/Secrier_Data")
load("SNV_Indel_Combined_COSMIC.Rda")
load("datOAC.clinicalPlusBFB.RData")

#merge the two
clinical_SNVIndel_data <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                                "HasBFB","BFBcateg")],
                                SNV_Indel_data[,c("sample","gene")],
                                by.x="Illumina_Barcode", by.y="sample",
                                all.x=TRUE, all.y=FALSE)

### Count the frequency of mutations across BFB types:
tb <- table(clinical_SNVIndel_data[,c("BFBtype","gene")])

#all BFB types
BFB_types <- sort(unique(clinical_SNVIndel_data$BFBtype))
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
BFBlike_data_SNVINDEL <- BFBtypes_df_list[1]
save(BFBlike_data_SNVINDEL, file="BFBlike_analysis_SNVINDEL.Rda")
DoubleMinute_data_SNVINDEL <- BFBtypes_df_list[2]
save(DoubleMinute_data_SNVINDEL, file="DoubleMinute_analysis_SNVINDEL.Rda")
FocalAmplifications_data_SNVINDEL <- BFBtypes_df_list[3]
save(FocalAmplifications_data_SNVINDEL, file="FocalAmplifications_analysis_SNVINDEL.Rda")
No_BFB_data_SNVINDEL <- BFBtypes_df_list[4]
save(No_BFB_data_SNVINDEL, file="NO_BFB_analysis_SNVINDEL.Rda")
Subtelomeric_data_SNVINDEL <- BFBtypes_df_list[5]
save(Subtelomeric_data_SNVINDEL, file = "Subtelomeric_analysis_SNVINDEL.Rda")
TandemDuplication_data_SNVINDEL <- BFBtypes_df_list[6]
save(TandemDuplication_data_SNVINDEL, file = "TandemDuplication_analysis_SNVINDEL.Rda")

### looking at results of above
###start here
setwd("~/Documents/R/Secrier_Data")
load("BFBlike_analysis_SNVINDEL.Rda")
load("DoubleMinute_analysis_SNVINDEL.Rda")
load("FocalAmplifications_analysis_SNVINDEL.Rda")
load("NO_BFB_analysis_SNVINDEL.Rda")
load("Subtelomeric_analysis_SNVINDEL.Rda")
load("TandemDuplication_analysis_SNVINDEL.Rda")

#BFB like
for (i in 1:length(BFBlike_data_SNVINDEL$`BFB-like`$pval_adjusted)){
  if(BFBlike_data_SNVINDEL$`BFB-like`$pval_adjusted[i] < 0.05){
    print(BFBlike_data_SNVINDEL$`BFB-like`$Gene_list[i])
    print(BFBlike_data_SNVINDEL$`BFB-like`$pval_adjusted[i])
  }
}
#CREBBP, p = 0.0007158703
#SGOL1, p = 0.0009638532

#Double Minute
for (i in 1:length(DoubleMinute_data_SNVINDEL$DoubleMinute$pval_adjusted)){
  if(DoubleMinute_data_SNVINDEL$DoubleMinute$pval_adjusted[i] < 0.05){
    print(DoubleMinute_data_SNVINDEL$DoubleMinute$Gene_list[i])
    print(DoubleMinute_data_SNVINDEL$DoubleMinute$pval_adjusted[i])
  }
}
#PRKAG1, p = 0.000002445238

#Focal Amplification
for (i in 1:length(FocalAmplifications_data_SNVINDEL$FocalAmplification$pval_adjusted)){
  if(FocalAmplifications_data_SNVINDEL$FocalAmplification$pval_adjusted[i] < 0.05){
    print(FocalAmplifications_data_SNVINDEL$FocalAmplification$Gene_list[i])
    print(FocalAmplifications_data_SNVINDEL$FocalAmplification$pval_adjusted[i])
  }
}
#TNNT3, 0.002630963

#None
for (i in 1:length(No_BFB_data_SNVINDEL$None$pval_adjusted)){
  if(No_BFB_data_SNVINDEL$None$pval_adjusted[i] < 0.05){
    print(No_BFB_data_SNVINDEL$None$Gene_list[i])
    print(No_BFB_data_SNVINDEL$None$pval_adjusted[i])
  }
}
#no significant results

#Subtelomeric
for (i in 1:length(Subtelomeric_data_SNVINDEL$Subtelomeric$pval_adjusted)){
  if(Subtelomeric_data_SNVINDEL$Subtelomeric$pval_adjusted[i] < 0.05){
    print(Subtelomeric_data_SNVINDEL$Subtelomeric$Gene_list[i])
    print(Subtelomeric_data_SNVINDEL$Subtelomeric$pval_adjusted[i])
  }
}
#CDKN2A, p = 0.00000009275232
#KIF20A, p = 0.002355915

#Tandem duplication
for (i in 1:length(TandemDuplication_data_SNVINDEL$TandemDuplication$pval_adjusted)){
  if(TandemDuplication_data_SNVINDEL$TandemDuplication$pval_adjusted[i] < 0.05){
    print(TandemDuplication_data_SNVINDEL$TandemDuplication$Gene_list[i])
    print(TandemDuplication_data_SNVINDEL$TandemDuplication$pval_adjusted[i])
  }
}
#no significant results
