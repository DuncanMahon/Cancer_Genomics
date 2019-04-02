### oncogene amplification linear model

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

#load in oncogene data
load("five_thresh_onco_amp.Rda")

#extract amplified oncogenes from cohort data
onco_fivep <- cohort_df[which(cohort_df$Gene_Name %in% five_thresh_onco_amp),]
onco_amp <- onco_fivep[which(onco_fivep$CNstatus=="Amp"),]
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
gene_list_onco <- unique(onco_fivep$Gene_Name)
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
#interpret this w/ Maria

# lets try stepAIC() for model refinement
library(MASS)
mod_step <- stepAIC(mod_oa,trace = FALSE)
mod_step$anova
#final model is different
summary(mod_step)
#MYC is significant here!

#side note:
#will try mod step model but without the intercept value
mod_step_v2 <- glm(BFBstatus ~ 0 + MYC + NFATC2 + CCND1 + CDK6 + TRRAP, data = training_oa, family = "binomial")
summary(mod_step_v2)
#AIC is much higher, similar pattern but no significance

# need to create testing_oa
testing_oa <- df_oa[-which(rownames(df_oa) %in% rownames(training_oa)),]
nrow(testing_oa)

## visualizing prediction results
results_prediction <- predict.glm(mod_step,testing_oa,type = "response")
str(results_prediction)
names(results_prediction)
#need to create BFB vector to match this
CR <- cohort_ID_BFB$HasBFB
df_oa <- cbind(df_oa,CR)
BFB_testing <- df_oa[which(rownames(df_oa) %in% names(results_prediction)),]
head(BFB_testing)
nrow(BFB_testing)
rownames(BFB_testing)
names(results_prediction)
BFB_testing$CR <- as.factor(BFB_testing$CR)
colnames(prediction_analysis) <- c("prediction","CR")
prediction_analysis <- as.data.frame(prediction_analysis)
prediction_analysis$CR <- as.factor(prediction_analysis$CR)
# create boxplot of probability distributions for BFBstatus
library(ggplot2)
plot <- ggplot(prediction_analysis, aes(x = CR, y = prediction,fill=CR)) +
  geom_boxplot()+
  scale_x_discrete(labels = c("No","Yes"))+
  stat_compare_means()+
  scale_fill_discrete(labels = c("No","Yes"))
setwd("~/Documents/R/Secrier_Plots")
ggsave(filename = "CR_prediction_result.pdf", plot = plot)

## count table for probability results
a <- 0
b <- 0
c <- 0
d <- 0
for (i in 1:nrow(prediction_analysis)){
  if (prediction_analysis$prediction[i] < 0.5 & prediction_analysis$BFB_status[i] == "0"){
    a <- a + 1
  } else if (prediction_analysis$prediction[i] < 0.5 & prediction_analysis$BFB_status[i] == "1") {
    b <- b + 1
  } else if (prediction_analysis$prediction[i] > 0.5 & prediction_analysis$BFB_status[i] == "0"){
    c <- c + 1
  } else {
    d <- d + 1
  }
}
a
b
c
d

prediction_results_matrix <- matrix(c(a,b,c,d),nrow = 2, ncol = 2,byrow = TRUE)
rownames(prediction_results_matrix) <- c("p < 0.5","p > 0.5")
colnames(prediction_results_matrix) <- c("noBFB","BFB")
prediction_results_matrix

