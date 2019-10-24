#  Gene Differentiation 
#  L.Dong
#  25-10-2019

# clear the environment
rm(list=ls())
# clear the console
cat('\014')


# load libraries
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(caret)
library(Rtsne)

# read in sample information, keeping those that are AD
samples <- read.csv('../filter-data/all.csv')
samples <- subset(samples, samples$clinical_group == 'AD')

# read in transcriptome data, change row names and remove the var 'X'
transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes time to load
rownames(transcriptome) <- transcriptome$X
transcriptome <- t(transcriptome[,-c(1)]) # samples need to be rows for Rtsne

tsne <- Rtsne(as.matrix(transcriptome), dims=2) # takes time
ggplot(as.data.frame(tsne$Y), aes(x=V1, y=V2, color=samples$lesional)) +
  geom_point(size=1.5)
