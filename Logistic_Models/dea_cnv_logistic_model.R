#logistic regression, CNV data
#from differential expression results
#try -ve log fold changes with amplifications with losses in CN

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

#load in differential expression results
load("de_gene_table3.Rda")
#select for negative log fold changes
#in this case they all are

#filter for our gene list
der_cnv <- cohort_df[which(cohort_df$Gene_Name %in% rownames(de_gene_table3)),]
#select for loss
der_cnv <- der_cnv[which(der_cnv$CNstatus=="Loss"),]
der_cnv_unique <- subset(der_cnv, select = c("Illumina_Barcode","HasBFB"))
der_cnv_unique <- unique(der_cnv_unique)

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
der_arr <- array(0, dim = c(length(cohort),length(unique(der_cnv$Gene_Name))))
t <- der_cnv
ID <- cohort
rownames(der_arr) <- ID
gene_list_der <- unique(der_cnv$Gene_Name)
colnames(der_arr) <- gene_list_der
for (i in rownames(der_arr)){
  for (j in colnames(der_arr)){
    if(nrow(t[which((t$Illumina_Barcode==i)&(t$Gene_Name==j)),])>0){
      der_arr[i,j] <- 1
    }
  }
}

#merge in BFB info next
df <- data.frame(der_arr)
df <- cbind(df,BFBstatus)

set.seed(46)
training_der <- df[sample(nrow(df),81),]
nrow(training_der)
#check for unique values
length(unique(rownames(training_der)))

###logistic regression on our data
mod_der <- glm(BFBstatus~.,data = training_der,family = "binomial")
summary(mod_der)

# lets try stepAIC()
library(MASS)
mod_step <- stepAIC(mod_der,trace = FALSE)
summary(mod_step)

#remove the intercept?
mod_der2 <- glm(BFBstatus~0+UGT2B15,data = training_der,family = "binomial")
summary(mod_der2)
#UGT2B15 significant alone

#try with TFF1 too?
mod_der3 <- glm(BFBstatus~0+UGT2B15+TFF1,data = training_der,family = "binomial")
summary(mod_der3)

#just TFF1?
mod_der4 <- glm(BFBstatus~0+TFF1,data = training_der,family = "binomial")
summary(mod_der4)
#nope

# need to create testing_oa
testing <- df[-which(rownames(df) %in% rownames(training_der)),]
nrow(testing)
predict.glm(mod_der2,testing,type = "response")
#this works
#these values are the probability of having these events

## visualizing prediction results
results_prediction <- predict.glm(mod_der2,testing,type = "response")
#need to create BFB vector to match this
CR <- cohort_ID_BFB$HasBFB
df <- cbind(df,CR)
BFB_testing <- df[which(rownames(df) %in% names(results_prediction)),]
head(BFB_testing)
nrow(BFB_testing)
rownames(BFB_testing)
names(results_prediction)
BFB_testing$CR <- as.factor(BFB_testing$CR)
prediction_analysis <- cbind(results_prediction,BFB_testing$CR)
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
ggsave(filename = "dif_exp_negfc_prediction_result.pdf", plot = plot)

## count table for probability results
a <- 0
b <- 0
c <- 0
d <- 0
for (i in 1:nrow(prediction_analysis)){
  if (prediction_analysis$prediction[i] < 0.5 & prediction_analysis$CR[i] == "0"){
    a <- a + 1
  } else if (prediction_analysis$prediction[i] < 0.5 & prediction_analysis$CR[i] == "1") {
    b <- b + 1
  } else if (prediction_analysis$prediction[i] > 0.5 & prediction_analysis$CR[i] == "0"){
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
