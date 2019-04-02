###finding associations between CN events and CR yes or no

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
                          cn_cg[,c("SampleID","CNstatus","LOH")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

###amplification
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

###deletion
yy2 <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$CNstatus == "Deletion")
yn2 <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$CNstatus != "Deletion")
ny2 <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$CNstatus == "Deletion")
nn2 <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$CNstatus != "Deletion")
n_yy2 <- length(yy2)
n_yn2 <- length(yn2)
n_ny2 <- length(ny2)
n_nn2 <- length(nn2)
#fishers exact test for deletion data
m2 <- rbind(c(n_yy2,n_ny2),c(n_yn2,n_nn2))
fisher.test(m2)

###LOH
yy3 <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$LOH == "yes")
yn3 <- which(dat.clinical.CNV$HasBFB == "Yes" & dat.clinical.CNV$LOH != "yes")
ny3 <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$LOH == "yes")
nn3 <- which(dat.clinical.CNV$HasBFB != "Yes" & dat.clinical.CNV$LOH != "yes")
n_yy3 <- length(yy3)
n_yn3 <- length(yn3)
n_ny3 <- length(ny3)
n_nn3 <- length(nn3)
#fishers exact test for LOH data
m3 <- rbind(c(n_yy3,n_ny3),c(n_yn3,n_nn3))
fisher.test(m3)

###deletion and loss combined
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

#creating a dataframe of our results
f <- fisher.test(m)
f2 <- fisher.test(m2)
f3 <- fisher.test(m3)
f4 <- fisher.test(m4)
f_pval <- f$p.value
f_pval2 <- f2$p.value
f_pval3 <- f3$p.value
f_pval4 <- f4$p.value
Genomic.Event <- c("Amplification","Deletion","Deletion.or.Loss","Loss.of.Heterozygosity")
P.values <- c(f_pval,f_pval2,f_pval4,f_pval3)
Adjusted.P.values <- p.adjust(P.values,method ="BH")
df_fisher <- data.frame(Genomic.Event,P.values,Adjusted.P.values)
#saving result
save(df_fisher, file = "CNV_BFB_Associations_Fishers.Rda")

### view result
setwd("~/Documents/R/Secrier_Data")
load("CNV_BFB_Associations_Fishers.Rda")
df_fisher

