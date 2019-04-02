### elastic net regularization
### Genes from MYC pathway (receptor tyrosine kinases)
### expression

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

#prepare data
m <- as.matrix(expr_rtks)
train_rows <- sample(1:36,36)
m_training <- m[train_rows,]
m_test <- m[-train_rows,]
BFB_train <- BFBstatus[train_rows]
BFB_test <- BFBstatus[-train_rows]

#load in package
#install.packages("glmnet")
library(glmnet)

#fit models
#??glmnet
fit.lasso <- glmnet(m_training, BFB_train, family="binomial", alpha=1)
fit.ridge <- glmnet(m_training, BFB_train, family="binomial", alpha=0)
fit.elnet <- glmnet(m_training, BFB_train, family="binomial", alpha=.5)

#20 fold cross validation
for(i in 0:20){
  assign(paste("fit",i,sep=""), cv.glmnet(m_training,BFB_train,type.measure = "mse",alpha=i/20,family="binomial"))
}

#plotting solutions
par(mfrow=c(3,2))
plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")

plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")

plot(fit.elnet, xvar="lambda")
plot(fit5, main="Elastic Net")

#MSE on test set
#?predict
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

#all values the same
#1.701438