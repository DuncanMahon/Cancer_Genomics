### add new columns for cnv file
## including one for presence/absence of loss of heterozygosity
## another for 'event' meaning 'amplification, neutral, deletion or loss'

#load in data
setwd("~/Documents/R/Secrier_Data")
load("copynumbers_129cohort.RData")
cn <- cndatatosave

#new column with CN status description
cn$CNstatus <- sapply(cn$AvTotalCN, 
                      function(x) ifelse(x>8,"Amp",
                                         ifelse(x>1,"Neutral",
                                                ifelse(x>0,"Loss","Deletion"))))

#new column w/ loss of heterozygosity
cn$LOH <- sapply(cn$AvMinorCN,
                 function(x) ifelse(x>0, "no",
                                    ifelse(x==0,"yes")))

#saving new data table
save(cn,file = "cn_129_cohort_updated_NOV52018")
