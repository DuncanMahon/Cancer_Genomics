###association between CNV and CR events over each gene

##standard start, load/merge
#Load my CN data
setwd("~/Documents/R/Secrier_Data")
load("cn_129_cohort_updated_NOV52018")

#reading in gene list
cancer_genes_df <- read.csv("COSMICCancerGeneCensus_allMon_Oct_8_13_35_11 2018.csv")
cancer_genes <- cancer_genes_df[,1]

#filtering cn df by cancer genes
cn_cg <- cn[which(cn$Gene_Name %in% cancer_genes),]

#reading clinical data
load("datOAC.clinicalPlusBFB.RData")

#merging clinical and CNV data
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn_cg[,c("SampleID","CNstatus","LOH","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

##now by gene analysis

#divide into amp, deletion, deletion OR loss, LOH
cln_cnv_amp <- dat.clinical.CNV[which(grepl("Amp",dat.clinical.CNV$CNstatus)),]
cln_cnv_del <- dat.clinical.CNV[which(grepl("Deletion",dat.clinical.CNV$CNstatus)),]
cln_cnv_loss <- dat.clinical.CNV[which(grepl("Loss",dat.clinical.CNV$CNstatus)),]
cln_cnv_delloss <- rbind(cln_cnv_del,cln_cnv_loss)
cln_cnv_LOH <- dat.clinical.CNV[which(grepl("yes",dat.clinical.CNV$LOH)),]

##amplification first
# Count the frequency of mutations in BFB vs non BFB samples:
tb_y_n_amp <- table(cln_cnv_amp[,c("HasBFB","Gene_Name")])

#hasbfb vs nobfb
no_BFB_sum_amp <- rowSums(tb_y_n_amp)[1]
yes_BFB_sum_amp <- rowSums(tb_y_n_amp)[2]
p_vals_amp <- c()
for (i in 1:ncol(tb_y_n_amp)){
  k <- tb_y_n_amp[,i]
  k_no <- k[1]
  k_yes <- k[2]
  other_no <- no_BFB_sum_amp - k_no
  other_yes <- yes_BFB_sum_amp - k_yes
  m <- rbind(c(k_yes,k_no),c(other_yes,other_no))
  f <- fisher.test(m)
  f_pval <- f$p.value
  p_vals_amp <- append(p_vals_amp,f_pval)
}
p_vals_amp_adjusted <- p.adjust(p_vals_amp,method="BH")


##deletion
# Count the frequency of mutations in BFB vs non BFB samples:
tb_y_n_del <- table(cln_cnv_del[,c("HasBFB","Gene_Name")])

#hasbfb vs nobfb
no_BFB_sum_del <- rowSums(tb_y_n_del)[1]
yes_BFB_sum_del <- rowSums(tb_y_n_del)[2]
p_vals_del <- c()
for (i in 1:ncol(tb_y_n_del)){
  k <- tb_y_n_del[,i]
  k_no <- k[1]
  k_yes <- k[2]
  other_no <- no_BFB_sum_del - k_no
  other_yes <- yes_BFB_sum_del - k_yes
  m <- rbind(c(k_yes,k_no),c(other_yes,other_no))
  f <- fisher.test(m)
  f_pval <- f$p.value
  p_vals_del <- append(p_vals_del,f_pval)
}
p_vals_del_adjusted <- p.adjust(p_vals_del,method="BH")


##deletion + loss
# Count the frequency of mutations in BFB vs non BFB samples:
tb_y_n_delloss <- table(cln_cnv_delloss[,c("HasBFB","Gene_Name")])

#hasbfb vs nobfb
no_BFB_sum_delloss <- rowSums(tb_y_n_delloss)[1]
yes_BFB_sum_delloss <- rowSums(tb_y_n_delloss)[2]
p_vals_delloss <- c()
for (i in 1:ncol(tb_y_n_delloss)){
  k <- tb_y_n_delloss[,i]
  k_no <- k[1]
  k_yes <- k[2]
  other_no <- no_BFB_sum_delloss - k_no
  other_yes <- yes_BFB_sum_delloss - k_yes
  m <- rbind(c(k_yes,k_no),c(other_yes,other_no))
  f <- fisher.test(m)
  f_pval <- f$p.value
  p_vals_delloss <- append(p_vals_delloss,f_pval)
}
p_vals_delloss_adjusted <- p.adjust(p_vals_delloss,method="BH")


##LOH
# Count the frequency of mutations in BFB vs non BFB samples:
tb_y_n_LOH <- table(cln_cnv_LOH[,c("HasBFB","Gene_Name")])

#hasbfb vs nobfb
no_BFB_sum_LOH <- rowSums(tb_y_n_LOH)[1]
yes_BFB_sum_LOH <- rowSums(tb_y_n_LOH)[2]
p_vals_LOH <- c()
for (i in 1:ncol(tb_y_n_LOH)){
  k <- tb_y_n_LOH[,i]
  k_no <- k[1]
  k_yes <- k[2]
  other_no <- no_BFB_sum_LOH - k_no
  other_yes <- yes_BFB_sum_LOH - k_yes
  m <- rbind(c(k_yes,k_no),c(other_yes,other_no))
  f <- fisher.test(m)
  f_pval <- f$p.value
  p_vals_LOH <- append(p_vals_LOH,f_pval)
}
p_vals_LOH_adjusted <- p.adjust(p_vals_LOH,method="BH")

##results, nothing to save, all adjusted p values are 1


