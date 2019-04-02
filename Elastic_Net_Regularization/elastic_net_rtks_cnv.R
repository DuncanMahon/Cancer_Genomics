### Elastic Net Regularization
### Genes from MYC pathway (receptor tyrosine kinases)
### Still looking at amplification of oncogenes
### CNV

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

#defining cohort data
cohort <- unique(dat.clinical.CNV$Illumina_Barcode)
cohort_df <- dat.clinical.CNV
#removing quotes
cohort$Gene_Name <- noquote(cohort_df$Gene_Name)

#load in rtk's
#install.packages("gdata")
require(gdata)
rtks <- read.xls("rtks.xlsx", sheet = 1, header = FALSE)

#extract amplified rtks from cohort data
onco_rtks <- cohort_df[which(noquote(cohort_df$Gene_Name) %in% rtks$V1),]
onco_amp <- onco_rtks[which(onco_rtks$CNstatus=="Amp"),]
onco_unique <- subset(onco_amp, select = c("Illumina_Barcode","HasBFB"))
onco_unique <- unique(onco_unique)

#convert Yes and No into 1 and 0 for BFBs
#selecting filtered ID list
BFBstatus <- c()
cohort_ID_BFB <- unique(subset(onco_amp,select = c("Illumina_Barcode","HasBFB")))
for (i in 1: length(cohort_ID_BFB$HasBFB)){
  if (cohort_ID_BFB$HasBFB[i] == "Yes"){
    BFBstatus <- append(BFBstatus, 1)
  } else {
    BFBstatus <- append(BFBstatus, 0)
  }
}

#Create and fill an array in preparation for linear modelling
#col dimension edited due to not having an amplificatino of any of these genes in some ID's
onco_arr <- array(0, dim = c(length(onco_amp$Illumina_Barcode),length(unique(onco_amp$Gene_Name))))
t <- onco_rtks
ID <- onco_amp$Illumina_Barcode
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

#leave it here for matrix form
#training rows now need to be edited
m <- as.matrix(onco_arr)
set.seed(46)
#49 is new 2/3rds
train_rows <- sample(1:49,49)
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

#plotting solutions (qc)
par(mfrow=c(3,2))
plot(fit.lasso, xvar="lambda")
plot(fit20, main="LASSO")

plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")

plot(fit.elnet, xvar="lambda")
plot(fit10, main="Elastic Net")

#MSE on test set
?predict
yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=m_test)
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=m_test)
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=m_test)
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=m_test)
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=m_test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=m_test)
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=m_test)
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=m_test)
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=m_test)
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=m_test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=m_test)
yhat11 <- predict(fit11, s=fit11$lambda.1se, newx=m_test)
yhat12 <- predict(fit12, s=fit12$lambda.1se, newx=m_test)
yhat13 <- predict(fit13, s=fit13$lambda.1se, newx=m_test)
yhat14 <- predict(fit14, s=fit14$lambda.1se, newx=m_test)
yhat15 <- predict(fit15, s=fit15$lambda.1se, newx=m_test)
yhat16 <- predict(fit16, s=fit16$lambda.1se, newx=m_test)
yhat17 <- predict(fit17, s=fit17$lambda.1se, newx=m_test)
yhat18 <- predict(fit18, s=fit18$lambda.1se, newx=m_test)
yhat19 <- predict(fit19, s=fit19$lambda.1se, newx=m_test)
yhat20 <- predict(fit20, s=fit20$lambda.1se, newx=m_test)

mse0 <- mean((BFB_test - yhat0)^2)
mse1 <- mean((BFB_test - yhat1)^2)
mse2 <- mean((BFB_test - yhat2)^2)
mse3 <- mean((BFB_test - yhat3)^2)
mse4 <- mean((BFB_test - yhat4)^2)
mse5 <- mean((BFB_test - yhat5)^2)
mse6 <- mean((BFB_test - yhat6)^2)
mse7 <- mean((BFB_test - yhat7)^2)
mse8 <- mean((BFB_test - yhat8)^2)
mse9 <- mean((BFB_test - yhat9)^2)
mse10 <- mean((BFB_test - yhat10)^2)
mse11 <- mean((BFB_test - yhat11)^2)
mse12 <- mean((BFB_test - yhat12)^2)
mse13 <- mean((BFB_test - yhat13)^2)
mse14 <- mean((BFB_test - yhat14)^2)
mse15 <- mean((BFB_test - yhat15)^2)
mse16 <- mean((BFB_test - yhat16)^2)
mse17 <- mean((BFB_test - yhat17)^2)
mse18 <- mean((BFB_test - yhat18)^2)
mse19 <- mean((BFB_test - yhat19)^2)
mse20 <- mean((BFB_test - yhat20)^2)

mse0
mse1
mse2
mse3
mse4
mse5
mse6
mse7
mse8
mse9
mse10
mse11
mse12
mse13
mse14
mse15
mse16
mse17
mse18
mse19
mse20

#no results given, ie NA




