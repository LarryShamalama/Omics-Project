rm(list=ls())
cat('\014')

library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(limma)
library(data.table)

samples <- read.csv('../filter-data/all.csv')
samples <- subset(samples, samples$clinical_group == 'AD')

transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes time to load
rownames(transcriptome) <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]

transcriptome <- transcriptome[rownames(transcriptome) %in% samples$sample_id,]
indices <- match(samples$sample_id, rownames(transcriptome))
stopifnot(sum(is.na(indices)) == 0)

transcriptome <- transcriptome[indices,]
stopifnot(all(rownames(transcriptome) == samples$sample_id))
transcriptome <- t(transcriptome)

lesion <- 2-as.integer(samples$lesional)
samples$lesional_new[samples$lesional=="LES"]<- 1
samples$lesional_new[samples$lesional=="NON_LES"]<- 0   
design <- as.data.frame(model.matrix(~lesion))

# fitting linear model
fit <- lmFit(transcriptome, design=design)
fit <- eBayes(fit, robust=TRUE)

signif <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif)