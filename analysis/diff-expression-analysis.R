#  Gene Differentiation 
#  L.Dong
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
#design<- as.data.frame(model.matrix(~lesional_new, data=samples ))
paired.design<- as.data.frame(model.matrix(~lesional_new + drop.levels(MAARS_identifier), data=samples))

# fitting linear model    NOTE:  getting a warning message about NAs
fit <- lmFit(transcriptome, design=paired.design)
# eBayes makes the variance more flexible
fit <- eBayes(fit, robust=TRUE)

# significance test
signif <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif)[ ,2]


#### ignore all this
# list of p-values for lesional/non-lesional per gene
gene_pvalues<- as.data.frame(fit$p.value[,2])
# calculate absolute value
gene_pvalues$absvalue<- abs(gene_pvalues$lesional_new)
# order by increasing abs value p-values
gene_pvalues<- order(gene_pvalues, by=absvalue)
####


#write genes to a csv file
write.fit(fit, signif, "myoutputdata.csv", adjust="BH")


