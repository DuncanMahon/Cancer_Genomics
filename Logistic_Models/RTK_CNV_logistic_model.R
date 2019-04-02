### logistic regression
### Genes from MYC pathway (receptor tyrosine kinases)
### Still looking at amplification of oncogenes
### CNV

## difference is we filter for preselected genes, instead of 
##greater than 5 threshold and we still filter for amplifications

### Load in dataset and merge
setwd("~/Documents/R/Secrier_Data")
load("datOAC.clinicalPlusBFB.RData")
load("cn_129_cohort_updated_NOV52018")
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn[,c("SampleID","CNstatus","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

#defining cohort data
cohort <- unique(dat.clinical.CNV$Illumina_Barcode)
cohort_df <- dat.clinical.CNV
#removing quotes
cohort_df$Gene_Name <- noquote(cohort_df$Gene_Name)

#load in rtk's
#install.packages("gdata")
require(gdata)
rtks <- read.xls("rtks.xlsx", sheet = 1, header = FALSE)

#extract amplified rtks from cohort data
onco_rtk <- cohort_df[which(cohort_df$Gene_Name %in% rtks$V1),]
onco_amp <- onco_rtk[which(onco_rtk$CNstatus=="Amp"),]
onco_unique <- subset(onco_amp, select = c("Illumina_Barcode","HasBFB"))
onco_unique <- unique(onco_unique)

#convert Yes and No into 1 and 0 for BFBs
BFBstatus <- c()
cohort_ID_BFB <- unique(subset(cohort_df,select = c("Illumina_Barcode","HasBFB")))
for (i in 1: length(cohort_ID_BFB$HasBFB)){
  if (cohort_ID_BFB$HasBFB[i] == "Yes"){
    BFBstatus <- append(BFBstatus, 1)
  } else {
    BFBstatus <- append(BFBstatus, 0)
  }
}

#Create and fill an array in preparation for linear modelling
onco_arr <- array(0, dim = c(length(cohort),length(unique(onco_amp$Gene_Name))))
t <- onco_amp
ID <- cohort
rownames(onco_arr) <- ID
gene_list_onco <- unique(onco_amp$Gene_Name)
colnames(onco_arr) <- gene_list_onco
for (i in rownames(onco_arr)){
  for (j in colnames(onco_arr)){
    if(nrow(t[which((t$Illumina_Barcode==i)&(t$Gene_Name==j)),])>0){
      onco_arr[i,j] <- 1
    }
  }
}

#merge in BFB info next
df_oa <- data.frame(onco_arr)
df_oa <- cbind(df_oa,BFBstatus)

#then take out 2/3 (randomly)
(nrow(df_oa)/3)*2
# will round to 81
set.seed(46)
training_oa <- df_oa[sample(nrow(df_oa),81),]
nrow(training_oa)
#check for unique values
length(unique(rownames(training_oa)))

###logistic regression on our data
mod_oa <- glm(BFBstatus~.,data = training_oa,family = "binomial")
summary(mod_oa)

# lets try stepAIC()
library(MASS)
mod_step <- stepAIC(mod_oa,trace = FALSE)
summary(mod_step)
#again MYC is significant

mod_oa2 <- glm(BFBstatus~MYC,data = training_oa,family = "binomial")
summary(mod_oa2)
#just MYC, not significant

mod_oa3 <- glm(BFBstatus~0+MYC,data = training_oa,family = "binomial")
summary(mod_oa3)
#minus the intercept, not significant
