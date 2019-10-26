#  Gene Differentiation 
#  L.Dong and K.Creel
#  24-10-2019

# clear the environment
rm(list=ls())
# clear the console
cat('\014')


# load libraries
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(rstudioapi)
library(limma)
library(data.table)

# read in sample information, keeping those that are AD
samples <- read.csv('../filter-data/all.csv')
samples <- subset(samples, samples$clinical_group == 'AD')

# read in transcriptome data, change row names and remove the var 'X'
transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes time to load
rownames(transcriptome) <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]

# keep the transcriptomes that are AD
transcriptome <- transcriptome[rownames(transcriptome) %in% samples$sample_id,]
# create an index, matching the sample_id to the row in the transcriptome
indices <- match(samples$sample_id, rownames(transcriptome))
# stop the process if there's any lines that have NA (sample_id has not matching transcriptome)
stopifnot(sum(is.na(indices)) == 0)

# order the transcriptome data according to the vector
transcriptome <- transcriptome[indices,]
# stop if not every row matches a sample_id
stopifnot(all(rownames(transcriptome) == samples$sample_id))
# transpose transcriptome data
transcriptome <- t(transcriptome)


#lesion <- 2-as.integer(samples$lesional)
samples$lesional_new[samples$lesional=="LES"]<- 1
samples$lesional_new[samples$lesional=="NON_LES"]<- 0   

library(gdata)
# create the model matrix for paired data
#   drop the unused dummy variables with drop.levels
paired.design<- as.data.frame(model.matrix(~lesional_new + drop.levels(MAARS_identifier), data=samples))

# fitting linear model
fit <- lmFit(transcriptome, design=paired.design)
# squeezes genewise residual variance towards a common value
fit <- eBayes(fit, robust=TRUE)

# significance tests - identify which genes are differentially expressed
signif <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif)[ ,2]


# extract top-ranked genes from the linear fit model
top.table <- topTable(fit, adjust.method="BH", p.value = 0.05, lfc = 1, sort.by = "P", n=Inf, coef=2)
gene_pvalues<- top.table[ ,(ncol(top.table)-3):ncol(top.table)]

# check to see how many genes are p<.05
length(which(top.table$adj.P.Val < 0.05))

# save gene pvalues
save(gene_pvalues,file="../filter-data/gene_pvalues.Rda")


# remove the _at at the end of gene (only do this when exporting, otherwise the gene name won't match the other files)
gene_names<- substr(row.names(gene_pvalues),1,nchar(row.names(gene_pvalues))-3)

# save gene names to a .csv
write.csv(rownames(gene_pvalues), "/Users/thewooz/Documents/Omics-Project/filter-data/genelist.csv")

