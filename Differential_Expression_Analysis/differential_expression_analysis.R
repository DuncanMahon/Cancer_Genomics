### expression analysis

#load in data
setwd("~/Documents/R/Secrier_Data")
load("mat.rawcounts.OAC.RData")
cts <- mat.rawcounts.OAC
#add BFB column
load("datOAC.clinicalPlusBFB.RData")
coldata <- cbind(dat.clinical$Illumina_Barcode,dat.clinical$HasBFB)
coldata <- as.data.frame(coldata)
colnames(coldata) <- c("Illumina_Barcode","HasBFB")

#filter coldata to only include IDs from expression analysis
coldata <- coldata[which(coldata$Illumina_Barcode %in% rownames(cts)),]
coldata <- unique(coldata)
#cross-check (should be 54), it is
length(coldata$Illumina_Barcode)

#transpose cts
cts <- t(cts)
#round each value
cts <- round(cts)

#DESeq2 dataset
library("DESeq2")
#create summarized experiment
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ HasBFB)
dds <- DESeq(dds)

#gain results
res <- results(dds)

#extracting gene lists for each p value
rownames(res)[which(res$padj < 0.05)]
rownames(res)[which(res$padj < 0.01)]
rownames(res)[which(res$padj < 0.001)]
#plenty at 0.05 to stay there

#create new gene table
#take top 50 absolute log fold changes from p < 0.05
genes2 <- rownames(res)[which(res$padj < 0.05)]
log2fc2 <- res$log2FoldChange[which(res$padj < 0.05)]
pval2 <- res$pvalue[which(res$padj < 0.05)]
padj2 <- res$padj[which(res$padj < 0.05)]
de_gene_table2 <- cbind(log2fc2,pval2,padj2)
rownames(de_gene_table2) <- genes2
top50 <- head((de_gene_table2[order(abs(log2fc2),decreasing = TRUE),]),n=50)

#create a heatmap
library("pheatmap")
#organise and normalize data
df <- as.data.frame(colData(dds)[,c("HasBFB")])
rownames(df) <- colnames(dds)
ntd <- normTransform(dds)

#create a for loop to create a pheatmamp going uo in intervals of ten so each instance
#can be visualised for clustering quality
s <- seq(10,150,10)
for (i in 1:length(s)){
  k <- head((de_gene_table2[order(abs(log2fc2),decreasing = TRUE),]),n=s[i])
  assign(paste("phm",s[i],sep = ""), 
         pheatmap(assay(ntd)[rownames(k),], 
                  cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df))
}

#stay at 50
phmap50 <- pheatmap(assay(ntd)[rownames(top50),], cluster_rows=TRUE, show_rownames=TRUE,
                    cluster_cols=TRUE, annotation_col=df, show_colnames = FALSE, 
                    annotation_legend = TRUE, annotation_names_col = FALSE)

#saving plot
setwd("~/Documents/R/Secrier_Plots")
pdf("phmap50.pdf")
phmap50
dev.off()

#saving our interesting cluster of genes
de_genes <- c("MUC6","TFF2","TFF1","ANXA10","NPC1L1","UGT2B15","SLC9A4")
setwd("~/Documents/R/Secrier_Data")
de_gene_table3 <- res[which(rownames(res) %in% de_genes),]
de_gene_table3 <- as.data.frame(de_gene_table3)
save(de_gene_table3,file="de_gene_table3.Rda")


