### elastic net regularization
### losses
### CNV

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

#load in del/loss data
load("five_thresh_deloss_tsg.Rda")

#extract lost TSGs from cohort data
tsg_fivep <- cohort_df[which(cohort_df$Gene_Name %in% five_thresh_deloss_tsg),]
tsg_del <- tsg_fivep[which(tsg_fivep$CNstatus=="Deletion"),]
tsg_unique <- subset(tsg_del, select = c("Illumina_Barcode","HasBFB"))
tsg_unique <- unique(tsg_unique)

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
tsg_arr <- array(0, dim = c(length(cohort),length(unique(tsg_del$Gene_Name))))
t <- tsg_del
ID <- cohort
rownames(tsg_arr) <- ID
gene_list_tsg <- unique(tsg_del$Gene_Name)
colnames(tsg_arr) <- gene_list_tsg
for (i in rownames(tsg_arr)){
  for (j in colnames(tsg_arr)){
    if(nrow(t[which((t$Illumina_Barcode==i)&(t$Gene_Name==j)),])>0){
      tsg_arr[i,j] <- 1
    }
  }
}

#leave it here for matrix form?
m <- as.matrix(tsg_arr)
set.seed(46)
train_rows <- sample(1:81,81)
m_training <- m[train_rows,]
m_test <- m[-train_rows,]
BFB_train <- BFBstatus[train_rows]
BFB_test <- BFBstatus[-train_rows]

#load in package
#install.packages("glmnet")
library(glmnet)

#fit models
??glmnet
fit.lasso <- glmnet(m_training, BFB_train, family="binomial", alpha=1)
fit.ridge <- glmnet(m_training, BFB_train, family="binomial", alpha=0)
fit.elnet <- glmnet(m_training, BFB_train, family="binomial", alpha=.5)

#20 fold cross validation
for(i in 0:20){
  assign(paste("fit",i,sep=""), cv.glmnet(m_training,BFB_train,type.measure = "mse",alpha=i/20,family="binomial"))
}
#won't work
#not enough deletion cases?

