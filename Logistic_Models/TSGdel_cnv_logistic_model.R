### TSG deletion linear model for CNV data

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

#load in TSG data
load("five_thresh_deloss_tsg.Rda")

#extract lost TSG from cohort data
TSG_fivep <- cohort_df[which(cohort_df$Gene_Name %in% five_thresh_deloss_tsg),]
TSG_del <- TSG_fivep[which(TSG_fivep$CNstatus=="Deletion"),]
#take unique instances
TSG_unique <- subset(TSG_del, select = c("Illumina_Barcode","HasBFB"))
TSG_unique <- unique(TSG_unique)

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
TSG_arr <- array(0, dim = c(length(cohort),length(unique(TSG_del$Gene_Name))))
t <- TSG_del
ID <- cohort
rownames(TSG_arr) <- ID
gene_list_TSG <- unique(TSG_del$Gene_Name)
colnames(TSG_arr) <- gene_list_TSG
for (i in rownames(TSG_arr)){
  for (j in colnames(TSG_arr)){
    if(nrow(t[which((t$Illumina_Barcode==i)&(t$Gene_Name==j)),])>0){
      TSG_arr[i,j] <- 1
    }
  }
}

#merge in BFB info next
df_TSG <- data.frame(TSG_arr)
df_TSG <- cbind(df_TSG,BFBstatus)

#then take out 2/3 (randomly)
(nrow(df_TSG)/3)*2
# will round to 81
set.seed(46)
training_TSG <- df_TSG[sample(nrow(df_TSG),81),]
nrow(training_TSG)
#check for unique values
length(unique(rownames(training_TSG)))

###logistic regression on our data
mod_TSG <- glm(BFBstatus~.,data = training_TSG,family = "binomial")
summary(mod_TSG)
#interpret this w/ Maria
#NAs are from where by chance the HasBFB examples have been removed

#stepAIC
library(MASS)
mod_step <- stepAIC(mod_TSG,trace = FALSE)
mod_step$anova
summary(mod_step)
#nothing significant


