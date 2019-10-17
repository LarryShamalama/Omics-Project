

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

transcriptome<- read.csv("../filter-data/transcriptome.csv")




str(features_all_samples)
str(gene_samples)


# Use the match() function to reorder the columns of the raw counts
reorder_idx <- match(rownames(features_all_samples), colnames(gene_samples))

# Reorder the columns of the count data
reordered_features <- gene_samples[ , reorder_idx_rmna]





library(edgeR)
library(limma)


# Model Matrix
design<- as.data.frame(model.matrix(~lesional, data=features_all_samples))
design

# Estimate common dispersion
gene_samples_commondisp<- edgeR::estimateCommonDisp(gene_samples)
gene_samples_commondisp

# Identify differentially expressed genes
fit <- glmFit(gene_samples_commondisp, design)
lrt <- glmLRT(fit, coef = 2)
signif_edgeR <- decideTestsDGE(lrt, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif_edgeR)




