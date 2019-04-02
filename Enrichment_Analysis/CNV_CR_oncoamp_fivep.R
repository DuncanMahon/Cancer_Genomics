### Oncogene/amp 5% and CR's

#load data
setwd("~/Documents/R/Secrier_Data")
load("cn_129_cohort_updated_NOV52018")
load("datOAC.clinicalPlusBFB.RData")
load("five_thresh_onco_amp.Rda")

#filtering cn df by 5% genes
cn_cog_five <- cn[which(cn$Gene_Name %in% five_thresh_onco_amp),]

#merging clinical and CNV data
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn_cog_five[,c("SampleID","CNstatus","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

#amplification of oncogenes
#extracting each element for fishers exact test
#in names, elemnt 1 refers to BFB or no BFB, elemtn 2 to Amp or no Amp
yy <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$CNstatus == "Amp")
yn <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$CNstatus != "Amp")
ny <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$CNstatus == "Amp")
nn <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$CNstatus != "Amp")

n_yy <- length(yy)
n_yn <- length(yn)
n_ny <- length(ny)
n_nn <- length(nn)

#fishers exact test for amplification data
m <- rbind(c(n_yy,n_ny),c(n_yn,n_nn))
fisher.test(m)
#result is 0.01174
#p < 0.05

###by each gene

#taking amplifications
cln_cnv_amp <- dat.clinical.CNV[which(grepl("Amp",dat.clinical.CNV$CNstatus)),]

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
# not quite all 1's but still not there
# each individual one is still over 0.05
