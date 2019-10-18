# Omics Presentation - Differential Expression
#  K. Creel
#  17 Oct 2019

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read in the cleaned data
transcriptome<- read.csv("../filter-data/transcriptome.csv")
description<- read.csv("../filter-data/all.csv")

library(edgeR)
library(limma)


## DIFFERENTIAL EXPRESSION

# in the 'description' file, create a new var for 'lesional' such that 0=non-lesion and 1=lesion
description$lesional_new[description$lesional=="LES"]<- 1
description$lesional_new[description$lesional=="NON_LES"]<- 0                      
                         
# create DGE object
AD_DGE<- DGEList(counts=transcriptome[,-1], 
                 samples=description, 
                 genes=transcriptome[,1],)

# Estimate common dispersion, add to the DGE object
AD_DGE_commondisp<- edgeR::estimateCommonDisp(AD_DGE)

# Model Matrix
design<- as.data.frame(model.matrix(~lesional_new, data=description))
head(design)

# Identify differentially expressed genes
fit <- glmFit(AD_DGE_commondisp, design)
lrt <- glmLRT(fit, coef = 2)
signif_edgeR <- decideTestsDGE(lrt, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif_edgeR)



## Using VOOM

# estimate the heteroscedasticity weights
v <- voom(transcriptome[,-1], design, plot = TRUE)

# fit the gaussian linear model with precision weights v
fit <- lmFit(v, design)

# Empirical Bayes test
fit <- eBayes(fit, robust = TRUE)
voomlimma_signif <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(voomlimma_signif)



