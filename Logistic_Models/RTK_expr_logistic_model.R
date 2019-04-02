### logistic regression
### Genes from MYC pathway (receptor tyrosine kinases)
### Still looking at amplification of oncogenes
### expression data

## difference is we filter for preselected genes, instead of 
##greater than 5 threshold and we don't filter for amplifications

### Load in dataset and merge
setwd("~/Documents/R/Secrier_Data")
load("datOAC.clinicalPlusBFB.RData")
load("cn_129_cohort_updated_NOV52018")
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn[,c("SampleID","CNstatus","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

### load in expression data
load("mat.expr.OAC.RData")
expr <- mat.expr.OAC

#load in rtk's
#install.packages("gdata")
require(gdata)
rtks <- read.xls("rtks.xlsx", sheet = 1, header = FALSE)

#create rtks df
expr_rtks <- expr[,which(colnames(expr) %in% rtks$V1)]

#extract and merge BFB data to expression data
cohort_df <- dat.clinical.CNV
cohort_df <- cohort_df[which(cohort_df$Illumina_Barcode %in% rownames(expr)),]
length(unique(cohort_df$Illumina_Barcode))
BFBstatus <- c()
cohort_ID_BFB <- unique(subset(cohort_df,select = c("Illumina_Barcode","HasBFB")))
for (i in 1: length(cohort_ID_BFB$HasBFB)){
  if (cohort_ID_BFB$HasBFB[i] == "Yes"){
    BFBstatus <- append(BFBstatus, 1)
  } else {
    BFBstatus <- append(BFBstatus, 0)
  }
}
expr_rtks_BFB <- cbind(expr_rtks,BFBstatus)

#randomly take out 2/3
(nrow(expr_rtks_BFB)/3)*2
set.seed(46)
training_expr <- expr_rtks_BFB[sample(nrow(expr_rtks_BFB),36),]
training_expr_df <- as.data.frame(training_expr)

#logistic regression on our data
mod_expr <- glm(BFBstatus~.,data = training_expr_df,family = "binomial")
summary(mod_expr)

#refine model by stepAIC
library(MASS)
mod_expr_step <- stepAIC(mod_expr,trace = FALSE)
summary(mod_expr_step)
#nothing significant


