#logistic regression, expression data
#from differential expression results

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

#load in dif expr results
load("de_gene_table3.Rda")
#create df of these genes
expr_der <- expr[,which(colnames(expr) %in% rownames(de_gene_table3))]

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
expr_der_BFB <- cbind(expr_der,BFBstatus)

#randomly take out 2/3
(nrow(expr_der_BFB)/3)*2
set.seed(46)
training_expr <- expr_der_BFB[sample(nrow(expr_der_BFB),36),]
training_expr_df <- as.data.frame(training_expr)

#logistic regression on our data
mod_expr <- glm(BFBstatus~.,data = training_expr_df,family = "binomial")
summary(mod_expr)

#refine model by stepAIC
library(MASS)
mod_expr_step <- stepAIC(mod_expr,trace = FALSE)
summary(mod_expr_step)
#nothing signficant in either model

#lets create expression boxplots for each gene
head(expr_der_BFB)
expr_der_BFB<-as.data.frame(expr_der_BFB)
expr_der_BFB$BFBstatus <- as.factor(expr_der_BFB$BFBstatus)
library("ggplot2")
library("ggpubr")
rownames(de_gene_table3)

p_ANXA10 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=ANXA10,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_MUC6 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=MUC6,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_NPC1L1 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=NPC1L1,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_SLC9A4 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=SLC9A4,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_TFF1 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=TFF1,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_TFF2 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=TFF2,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

p_UGT2B15 <- ggplot(expr_der_BFB,aes(x=BFBstatus,y=UGT2B15,group=BFBstatus,fill=BFBstatus))+
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means() +
  scale_fill_discrete(name="CR",
                      breaks = c("0","1"),
                      labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("0","1"),
                   labels = c("No","Yes")) +
  xlab("CR")

# saving all plots
setwd("~/Documents/R/Secrier_Plots")
ggsave(filename = "expr_ANXA10.jpg", plot = p_ANXA10)
ggsave(filename = "expr_MUC6.jpg", plot = p_MUC6)
ggsave(filename = "expr_NPC1L1.jpg", plot = p_NPC1L1)
ggsave(filename = "expr_SLC9A4.jpg", plot = p_SLC9A4)
ggsave(filename = "expr_TFF1.jpg", plot = p_TFF1)
ggsave(filename = "expr_TFF2.jpg", plot = p_TFF2)
ggsave(filename = "expr_UGT2B15.jpg", plot = p_UGT2B15)

#combine these plots
library(cowplot)
expression_plot <- plot_grid(p_ANXA10,
                             p_MUC6,
                             p_NPC1L1,
                             p_SLC9A4,
                             p_TFF1,
                             p_TFF2,
                             p_UGT2B15,
                             labels = c("A","B","C","D","E","F","G"),
                             nrow = 4,
                             ncol = 2)
ggsave(filename = "expression_plot2.jpg", plot = expression_plot)




