### survival analysis
### kaplan meier curve
### using our clinical data

#load in packages
library("survminer")
require("survival")

# Load clinical data
setwd("~/Documents/R/Secrier_Data")
load("datOAC.clinicalPlusBFB.RData")

##Fit survival curves
#BFB yes or no
colnames(dat.clinical)[which(names(dat.clinical) == "HasBFB")] <- "HasCR"
fit<- survfit(Surv(OS.Months, Cnsr) ~ HasCR, data = dat.clinical)
surv_p_HasCR <- ggsurvplot(fit, 
                            pval = TRUE, conf.int = TRUE,
                            risk.table = TRUE,
                            risk.table.col = "strata",
                            linetype = "strata",
                            ggtheme = theme_bw(),
                            palette = c("#4897D8", "#FFDB5C"))

#BFB category
colnames(dat.clinical)[which(names(dat.clinical) == "BFBcateg")] <- "CRcategory"
fit2<- survfit(Surv(OS.Months, Cnsr) ~ CRcategory, data = dat.clinical)
surv_p_CRcategory <- ggsurvplot(fit2, 
                                pval = TRUE, conf.int = TRUE,
                                risk.table = TRUE,
                                risk.table.col = "strata",
                                linetype = "strata",
                                ggtheme = theme_bw(),
                                palette = c("#4897D8", "#FFDB5C", "#FA6E59"))

#BFB type
fit3<- survfit(Surv(OS.Months, Cnsr) ~ BFBtype, data = dat.clinical)
ggsurvplot(fit3, 
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           ggtheme = theme_bw(),
           palette = c("#4897D8", "#FFDB5C", "#FA6E59","#F8A055","#E7B800", "#2E9FDF"))
#error, potentially does not like so many factors

#focus on focal amplifications
#creating a column of focal amplifications, yes or no
CR_FA <- rep("FocalAmplification",129)
for (i in 1:nrow(dat.clinical)) {
  if(dat.clinical$BFBtype[i] != "FocalAmplification"){
    CR_FA[i] <- "Other"
  }
}
table(CR_FA)
table(dat.clinical$BFBtype)

#focal amplifications
dat.clinical <- cbind(dat.clinical,CR_FA)
fit4<- survfit(Surv(OS.Months, Cnsr) ~ CR_FA, data = dat.clinical)
surv_p_CRfocalamplification <- ggsurvplot(fit4, 
                                          pval = TRUE, conf.int = TRUE,
                                          risk.table = TRUE,
                                          risk.table.col = "strata",
                                          linetype = "strata",
                                          ggtheme = theme_bw(),
                                          palette = c("#4897D8", "#FFDB5C"))


#focus on BFB-like
#creating a column of BFB-like, yes or no
BFB_like <- rep("BFB-like",129)
for (i in 1:nrow(dat.clinical)) {
  if(dat.clinical$BFBtype[i] != "BFB-like"){
    BFB_like[i] <- "Other"
  }
}
table(BFB_like)
table(dat.clinical$BFBtype)

#BFB-like
dat.clinical <- cbind(dat.clinical,BFB_like)
fit5<- survfit(Surv(OS.Months, Cnsr) ~ BFB_like, data = dat.clinical)
surv_p_BFBlike <- ggsurvplot(fit5, 
                             pval = TRUE, conf.int = TRUE,
                             risk.table = TRUE,
                             risk.table.col = "strata",
                             linetype = "strata",
                             ggtheme = theme_bw(),
                             palette = c("#4897D8", "#FFDB5C"))

#focus on tandem-duplication
#creating a column of tandem duplications, yes or no
CR_TD <- rep("TandemDuplication",129)
for (i in 1:nrow(dat.clinical)) {
  if(dat.clinical$BFBtype[i] != "TandemDuplication"){
    CR_TD[i] <- "Other"
  }
}
table(CR_TD)
table(dat.clinical$BFBtype)

#tandem duplication
dat.clinical <- cbind(dat.clinical,CR_TD)
fit6<- survfit(Surv(OS.Months, Cnsr) ~ CR_TD, data = dat.clinical)
surv_p_CRtandemduplication <- ggsurvplot(fit6, 
                                         pval = TRUE, conf.int = TRUE,
                                         risk.table = TRUE,
                                         risk.table.col = "strata",
                                         linetype = "strata",
                                         ggtheme = theme_bw(),
                                         palette = c("#4897D8", "#FFDB5C"))

#save these plots
setwd("~/Documents/R/Secrier_Plots")
surv_p_HasCR
surv_p_BFBlike
surv_p_CRcategory
surv_p_CRfocalamplification
surv_p_CRtandemduplication
#plots saved in GUI as for ggsurvplot
#ggsave() and jpeg() do not work for ggsurvplot


