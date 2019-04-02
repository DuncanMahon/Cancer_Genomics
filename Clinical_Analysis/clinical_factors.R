### Investigating clinical factors of CR status

### Load in dataset and merge
setwd("~/Documents/R/Secrier_Data")
load("datOAC.clinicalPlusBFB.RData")
#visualise clinical factors
colnames(dat.clinical)

#load relevant packages
library("ggplot2")
library("ggpubr")
library(cowplot)

# age
# (of diagnosis)
cd_age <- dat.clinical[,c("AGE.diagn","HasBFB")]
cd_age$AGE.diagn <- as.numeric(cd_age$AGE.diagn)
colnames(cd_age) <- c("Age.diagn","CR")
p_age <- ggplot(cd_age,aes(x=CR,y=Age.diagn,group=CR,fill=CR))+
  geom_boxplot() +
  scale_y_continuous(limits = c(40,90),breaks = seq(40,90,5)) +
  geom_jitter() +
  stat_compare_means() +
  ylab("Age of Diagnosis") +
  xlab("CR")

# drinking
# units
cd_units <- dat.clinical[,c("EX.Units.Per.Week.In.Total","HasBFB")]
cd_units$EX.Units.Per.Week.In.Total <- as.numeric(cd_units$EX.Units.Per.Week.In.Total)
colnames(cd_units) <- c("EX.Units.Per.Week.In.Total","CR")
p_units <- ggplot(cd_units,aes(x=CR,y=EX.Units.Per.Week.In.Total,fill=CR))+
  geom_boxplot() +
  scale_y_continuous(limits = c(0,80),breaks = seq(0,80,5)) +
  stat_compare_means() +
  geom_jitter() +
  ylab("Daily Units of Alcohol") +
  xlab("CR")


# smoking
# Average.Cigarrettes.Per.Day
cd_smoking <- dat.clinical[,c("Average.Cigarrettes.Per.Day","HasBFB")]
cd_smoking$Average.Cigarrettes.Per.Day <- as.numeric(cd_smoking$Average.Cigarrettes.Per.Day)
colnames(cd_smoking) <- c("Average.Cigarrettes.Per.Day","CR")
p_smoking <- ggplot(cd_smoking,aes(x=CR,y=Average.Cigarrettes.Per.Day,group=CR,fill=CR))+
  geom_boxplot() +
  scale_y_continuous(limits = c(0,60),breaks = seq(0,60,5)) +
  stat_compare_means() +
  geom_jitter() +
  ylab("Daily Cigarette Use")

# Barretts (pre EAC condition)
cd_barretts <- dat.clinical[,c("Barretts.y.n","HasBFB")]
table(cd_barretts)
cd_barretts <- cd_barretts[-which(cd_barretts$Barretts.y.n=="x"),]
table(cd_barretts)
fisher.test(table(cd_barretts))
# p = 0.3319

# sex
cd_sex <- dat.clinical[,c("Sex","HasBFB")]
table(cd_sex)
fisher.test(table(cd_sex))
# p = 0.8066

# S17A
# mutational signature
# these are continuous!
cd_s17A <- dat.clinical[,c("S17A","HasBFB")]
colnames(cd_s17A) <- c("S17A","CR")
p_s17A <- ggplot(cd_s17A,aes(x=CR,y=S17A,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# S2 (APOBEC)
#needs renaming
cd_s2 <- dat.clinical[,c("S2 (APOBEC)","HasBFB")]
colnames(cd_s2) <- c("S2","CR")
p_s2 <- ggplot(cd_s2,aes(x=CR,y= S2,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# S1 (age)
cd_s1 <- dat.clinical[,c("S1 (age)","HasBFB")]
colnames(cd_s1) <- c("S1","CR")
p_s1 <- ggplot(cd_s1,aes(x=CR,y= S1,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# S18-like
cd_s18 <- dat.clinical[,c("S18-like","HasBFB")]
colnames(cd_s18) <- c("S18","CR")
p_s18 <- ggplot(cd_s18,aes(x=CR,y=S18,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# S3 (BRCA)
cd_s3 <- dat.clinical[,c("S3 (BRCA)","HasBFB")]
colnames(cd_s3) <- c("S3","CR")
p_s3 <- ggplot(cd_s3,aes(x=CR,y=S3,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# S17B
cd_s17B <- dat.clinical[,c("S17B","HasBFB")]
colnames(cd_s17B) <- c("S17B","CR")
p_s17B <- ggplot(cd_s17B,aes(x=CR,y=S17B,group=CR,fill=CR))+
  geom_boxplot() +
  stat_compare_means() +
  geom_jitter()

# saving all plots
setwd("~/Documents/R/Secrier_Plots")
ggsave(filename = "BFB_age.jpg", plot = p_age)
ggsave(filename = "BFB_smoking.jpg", plot = p_smoking)
ggsave(filename = "BFB_units.jpg", plot = p_units)
ggsave(filename = "BFB_s1.jpg", plot = p_s1)
ggsave(filename = "BFB_s17A.jpg", plot = p_s17A)
ggsave(filename = "BFB_s17B.jpg", plot = p_s17B)
ggsave(filename = "BFB_s18.jpg", plot = p_s18)
ggsave(filename = "BFB_s2.jpg", plot = p_s2)
ggsave(filename = "BFB_s3.jpg", plot = p_s3)

#combining plots
clinical_plot <- plot_grid(p_age,
          p_smoking,
          p_units,
          p_s1,
          p_s2,
          p_s3,
          p_s17A,
          p_s17B,
          p_s18,
          labels = c('A','B','C','D','E','F','G','H','I'),
          label_x = 0.2,
          nrow = 3,
          ncol = 3)
#saving combined plot
ggsave(filename = "clinical_plot.jpg", plot = clinical_plot)




