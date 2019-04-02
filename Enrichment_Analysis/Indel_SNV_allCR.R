#Indel + SNV analysis (all BFB)

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

### Count the frequency of mutations in BFB vs non BFB samples:
tb_y_n <- table(clinical_SNVIndel_data[,c("HasBFB","gene")])
tb <- table(clinical_SNVIndel_data[,c("BFBtype","gene")])

#hasbfb vs nobfb
no_BFB_sum <- rowSums(tb_y_n)[1]
yes_BFB_sum <- rowSums(tb_y_n)[2]
p_vals <- c()
for (i in 1:ncol(tb_y_n)){
  k <- tb_y_n[,i]
  k_no <- k[1]
  k_yes <- k[2]
  other_no <- no_BFB_sum - k_no
  other_yes <- yes_BFB_sum - k_yes
  m <- rbind(c(k_yes,k_no),c(other_yes,other_no))
  f <- fisher.test(m)
  f_pval <- f$p.value
  p_vals <- append(p_vals,f_pval)
}
p_vals
pvals_adjusted <- p.adjust(p_vals,method="BH")

for (i in 1:length(pvals_adjusted)){
  if(pvals_adjusted[i] < 0.05){
    print(pvals_adjusted[i])
  }
}

###everything is 1 after adjustment