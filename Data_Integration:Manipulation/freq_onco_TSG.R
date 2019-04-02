### gene frequency summary
### frequency in alterations (a) mutations/indels, (b) amp, (c) del/loss
### in (1) tumour suppressors and (2) oncogenes

#load data
setwd("~/Documents/R/Secrier_Data")
load("COSMIC_TSG.Rda")
load("COSMIC_oncogenes.Rda")
load("cn_129_cohort_updated_NOV52018")
load("datOAC.clinicalPlusBFB.RData")
load("SNV_Indel_Combined_COSMIC.Rda")

#extract gene list
cancer_oncogenes <- can_gene_onco_all[,1]
cancer_TSG <- can_gene_TSG_all[,1]

#filtering cn df by cancer genes
cn_cog <- cn[which(cn$Gene_Name %in% cancer_oncogenes),]
#merging clinical and CNV data
dat.clinical.CNV <- merge(dat.clinical[,c("Illumina_Barcode","BFBtype",
                                          "HasBFB","BFBcateg")],
                          cn_cog[,c("SampleID","CNstatus","Gene_Name")],
                          by.x="Illumina_Barcode", by.y="SampleID",
                          all.x = TRUE, all.y = FALSE)

#taking amplifications
cln_cnv_amp <- dat.clinical.CNV[which(grepl("Amp",dat.clinical.CNV$CNstatus)),]
clin_cnv_amp_onco <- cln_cnv_amp[which(cln_cnv_amp$Gene_Name %in% cancer_oncogenes),]
clin_cnv_amp_tsg <- cln_cnv_amp[which(cln_cnv_amp$Gene_Name %in% cancer_TSG),]

#creating tables to summarise the presence of amplifications of oncogenes and
# deletion/loss of TSG's
#onco tables
tb_clin_cnv_amp_onco <- table(clin_cnv_amp_onco[,c("Gene_Name")])
tb_clin_cnv_amp_onco_BFB_nBFB <- table(clin_cnv_amp_onco[,c("HasBFB","Gene_Name")])

#TSG tables
tb_clin_cnv_deloss_tsg <- table(clin_cnv_deloss_tsg[,c("Gene_Name")])
tb_clin_cnv_deloss_tsg_BFB_nBFB <- table(clin_cnv_deloss_tsg[,c("HasBFB","Gene_Name")])

library(ggplot2)

##create dataframes for plotting
##also order them by frequency of amplified oncogenes and deletion/loss of TSG's respectively
cols <- c("Gene","Frequency")
df_clin_cnv_amp_onco <- as.data.frame(tb_clin_cnv_amp_onco)
colnames(df_clin_cnv_amp_onco) <- cols
df_clin_cnv_amp_onco <- df_clin_cnv_amp_onco[order(df_clin_cnv_amp_onco$Frequency, decreasing = TRUE),]
df_clin_cnv_amp_onco$Gene <- factor(df_clin_cnv_amp_onco$Gene,levels = df_clin_cnv_amp_onco$Gene[order(-df_clin_cnv_amp_onco$Frequency)])

df_clin_cnv_deloss_tsg <- as.data.frame(tb_clin_cnv_deloss_tsg)
colnames(df_clin_cnv_deloss_tsg) <- cols
df_clin_cnv_deloss_tsg <- df_clin_cnv_deloss_tsg[order(df_clin_cnv_deloss_tsg$Frequency, decreasing = TRUE),]
df_clin_cnv_deloss_tsg$Gene <- factor(df_clin_cnv_deloss_tsg$Gene,levels = df_clin_cnv_deloss_tsg$Gene[order(-df_clin_cnv_deloss_tsg$Frequency)])

### lets filter off and take oncogene amplifications over 5% and tsg del/loss over 5% of cohort
# first we want to take our critical value that for over 5% of our cohort (129 people are in our cohort)
# will take the figure for BFB and nBFB combined
k <- ceiling(129/20)

#filter for onco/amp
five_thresh_onco_amp <- df_clin_cnv_amp_onco$Gene[which(df_clin_cnv_amp_onco$Frequency >= k)]
save(five_thresh_onco_amp,file = "five_thresh_onco_amp.Rda")

#filter for tsg/deloss
five_thresh_deloss_tsg <- df_clin_cnv_deloss_tsg$Gene[which(df_clin_cnv_deloss_tsg$Frequency >= k)]
save(five_thresh_deloss_tsg,file = "five_thresh_deloss_tsg.Rda")

#lets also filter out top 10% for more stringent analysis
j <- ceiling(129/10)

#filter for onco/amp
ten_thresh_onco_amp <- df_clin_cnv_amp_onco$Gene[which(df_clin_cnv_amp_onco$Frequency >= j)]
save(ten_thresh_onco_amp,file = "ten_thresh_onco_amp.Rda")

#filter for tsg/deloss
ten_thresh_deloss_tsg <- df_clin_cnv_deloss_tsg$Gene[which(df_clin_cnv_deloss_tsg$Frequency >= j)]
#save(ten_thresh_deloss_tsg,file = "ten_thresh_deloss_tsg.Rda")

#create a plot to visualise most frequently amplified oncogenes and deleted/lost TSG's

#top 5% and 10% of amplified oncogenes and deleted/lost TSG's
cutoff_5 <- data.frame( x = c(-Inf, Inf), y = k, cutoff_5 = factor(k) )
cutoff_10 <- data.frame( x = c(-Inf, Inf), y = j, cutoff_10 = factor(j) )

onco_amp_5_10 <- ggplot(df_clin_cnv_amp_onco[1:20,], aes(x=Gene,y=Frequency)) +
  geom_bar(stat = "identity",fill="blue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_line(aes( x, y, linetype = cutoff_5 ), cutoff_5) +
  geom_line(aes( x, y, linetype = cutoff_10 ), cutoff_10) +
  scale_linetype(name = "Cutoffs",
                 breaks = c("13","7"),
                 labels=c("10%","5%"))
  

tsg_deloss_5_10 <- ggplot(df_clin_cnv_deloss_tsg[1:20,], aes(x=Gene,y=Frequency)) +
  geom_bar(stat = "identity",fill="orange") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_line(aes( x, y, linetype = cutoff_5 ), cutoff_5) +
  geom_line(aes( x, y, linetype = cutoff_10 ), cutoff_10) +
  scale_linetype(name = "Cutoffs",
                 breaks = c("13","7"),
                 labels=c("10%","5%"))

#combine the two plots
library(cowplot)
p_gene_summary <- plot_grid(onco_amp_5_10,
                            tsg_deloss_5_10,
                            labels = c("A","B"),
                            nrow = 2,
                            ncol = 1)
setwd("~/Documents/R/Secrier_Plots")
ggsave(filename = "gene_summary.jpg", plot = p_gene_summary)



