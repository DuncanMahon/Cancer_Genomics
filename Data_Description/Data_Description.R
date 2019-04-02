### Basic data description for lab report

### Load in dataset and merge
setwd("~/Documents/R/Secrier_Data")
load("datOAC.clinicalPlusBFB.RData")
load("cn_129_cohort_updated_NOV52018")

#visualise data
head(cn)
head(dat.clinical)

#find length of unique list of ID's
length(unique(dat.clinical$Illumina_Barcode))
#122 individuals have clinical description
length(unique(cn$Illumina_Barcode))
#129 in this case
#note they do not fully match (upstream error)

unique(dat.clinical$BFBtype)
#None; Focal Amplification; BFB-like; Tandem Duplication; Double Minute; Subtelomeric  
unique(dat.clinical$BFBcateg)
#None; Amplified; BFBlike

table(dat.clinical$Sex)
#21 female, 108 male
table(dat.clinical$STAGE)
#1-4, 36, 18, 43 and 22 respectively
#9 unknown
table(dat.clinical$Barretts.y.n)
#45 no and 72 yes respectively
#12 unknown
