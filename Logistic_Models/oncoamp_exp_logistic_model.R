### Logistic model using expression data
### oncogenic amplifications only

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

#load in oncogene data
load("ten_thresh_onco_amp.Rda")
#create onco df
expr_onco <- expr[,which(colnames(expr) %in% ten_thresh_onco_amp)]

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
expr_onco_BFB <- cbind(expr_onco,BFBstatus)

#randomly take out 2/3
(nrow(expr_onco_BFB)/3)*2
set.seed(46)
training_expr <- expr_onco_BFB[sample(nrow(expr_onco_BFB),36),]
training_expr_df <- as.data.frame(training_expr)

#logistic regression on our data
mod_expr <- glm(BFBstatus~.,data = training_expr_df,family = "binomial")
summary(mod_expr)

#refine model by stepAIC
library(MASS)
mod_expr_step <- stepAIC(mod_expr,trace = FALSE)
summary(mod_expr_step)
#nothing significant
#stepAIC only keeps ERBB2 in the model

#expression model, just MYC
mod_expr_MYC <- glm(BFBstatus~ 0 + MYC,data = training_expr_df,family="binomial")
summary(mod_expr_MYC)
#significant

#create boxplot for each gene, starting with MYC for expression
# between BFB 0 and BFB 1

#load packages
library("ggplot2")
library("ggpubr")

#create dataframe
expr_all_df <- as.data.frame(expr_onco_BFB)

#create each plot
cd_MYC <- cbind(expr_all_df[,c("MYC","BFBstatus")],cohort_ID_BFB$HasBFB)
colnames(cd_MYC) <- c("MYC","CR_binary","CR")
p_MYC <- ggplot(cd_MYC,aes(x=CR,y=MYC,group=CR,fill=CR))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

cd_CCND1 <- cbind(expr_all_df[,c("CCND1","BFBstatus")],cohort_ID_BFB$HasBFB)
colnames(cd_CCND1) <- c("CCND1","CR_binary","CR")
p_CCND1 <- ggplot(cd_CCND1,aes(x=CR,y=CCND1,group=CR,fill=CR))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

cd_CCR7 <- cbind(expr_all_df[,c("CCR7","BFBstatus")],cohort_ID_BFB$HasBFB)
colnames(cd_CCR7) <- c("CCR7","CR_binary","CR")
p_CCR7 <- ggplot(cd_CCR7,aes(x=CR,y=CCR7,group=CR,fill=CR))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

cd_CDK6 <- cbind(expr_all_df[,c("CDK6","BFBstatus")],cohort_ID_BFB$HasBFB)
colnames(cd_CDK6) <- c("CDK6","CR_binary","CR")
p_CDK6 <- ggplot(cd_CDK6,aes(x=CR,y=CDK6,group=CR,fill=CR))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

cd_ERBB2 <- cbind(expr_all_df[,c("ERBB2","BFBstatus")],cohort_ID_BFB$HasBFB)
colnames(cd_ERBB2) <- c("ERBB2","CR_binary","CR")
p_ERBB2 <- ggplot(cd_ERBB2,aes(x=CR,y=ERBB2,group=CR,fill=CR))+
  geom_boxplot()+
  geom_jitter()+
  stat_compare_means()

# saving all plots
setwd("~/Documents/R/Secrier_Plots")
ggsave(filename = "expr_MYC.jpg", plot = p_MYC)
ggsave(filename = "expr_CCND1.jpg", plot = p_CCND1)
ggsave(filename = "expr_CCR7.jpg", plot = p_CCR7)
ggsave(filename = "expr_CDK6.jpg", plot = p_CDK6)
ggsave(filename = "expr_ERBB2.jpg", plot = p_ERBB2)

#combine these plots
library(cowplot)
expression_plot <- plot_grid(p_MYC,
                             p_CCND1,
                             p_CCR7,
                             p_CDK6,
                             p_ERBB2,
                             labels = c("A","B","C","D","E"),
                             label_x = 0.2,
                             nrow = 2,
                             ncol = 3)
ggsave(filename = "expression_plot.jpg", plot = expression_plot)



