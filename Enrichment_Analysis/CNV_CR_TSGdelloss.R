###del/loss of TSG and CR

#load data
setwd("~/Documents/R/Secrier_Data")
load("COSMIC_TSG.Rda")
load("cn_129_cohort_updated_NOV52018")
load("datOAC.clinicalPlusBFB.RData")

#extract gene list
cancer_TSG <- can_gene_TSG_all[,1]

#filtering cn df by TSG
cn_tsg <- cn[which(cn$Gene_Name %in% cancer_TSG),]

#merging clinical and CNV data
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn_tsg[,c("SampleID","CNstatus","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

#deletion and loss combined
yy4 <- which(dat.clinical.CNV$HasBFB == "Yes" & (dat.clinical.CNV$CNstatus == "Deletion" | dat.clinical.CNV$CNstatus == "Loss"))
yn4 <- which(dat.clinical.CNV$HasBFB == "Yes" & (dat.clinical.CNV$CNstatus != "Deletion" & dat.clinical.CNV$CNstatus != "Loss"))
ny4 <- which(dat.clinical.CNV$HasBFB != "Yes" & (dat.clinical.CNV$CNstatus == "Deletion" | dat.clinical.CNV$CNstatus == "Loss"))
nn4 <- which(dat.clinical.CNV$HasBFB != "Yes" & (dat.clinical.CNV$CNstatus != "Deletion" & dat.clinical.CNV$CNstatus != "Loss"))
n_yy4 <- length(yy4)
n_yn4 <- length(yn4)
n_ny4 <- length(ny4)
n_nn4 <- length(nn4)
#fishers exact test for deletion data
m4 <- rbind(c(n_yy4,n_ny4),c(n_yn4,n_nn4))
fisher.test(m4)
#super significant, p < 0.0001

###now by gene

#filter deletions/losses
cln_cnv_del <- dat.clinical.CNV[which(grepl("Deletion",dat.clinical.CNV$CNstatus)),]
cln_cnv_loss <- dat.clinical.CNV[which(grepl("Loss",dat.clinical.CNV$CNstatus)),]
cln_cnv_delloss <- rbind(cln_cnv_del,cln_cnv_loss)

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
#same same, just = 1

