###adding sample ID to CNV file

#load in data
setwd("~/Documents/R/Secrier_Data")
load("cn_129_cohort_updated_NOV52018")

#split illumina barcode
cn$SampleID <- sapply(cn$Illumina_Barcode,
                      function(x) strsplit(x,"_vs_")[[1]][1])

#saving new data table
save(cn,file = "cn_129_cohort_updated_NOV52018")
