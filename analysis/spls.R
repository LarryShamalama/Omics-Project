# reset console
cat('\014')
rm(list=ls())

packages <- c('ggplot2',
              'rstudioapi',
              'data.table',
              'mixOmics',
              'SuperLearner',
              'parallel')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in packages){
  tryCatch(
    if (!(pkg %in% rownames(installed.packages()))){
      install.packages(pkg)
    }, error=function(e){
      BiocManager::install(pkg) # try installing with BiocManager
    })
}

library('ggplot2')
library('rstudioapi')
library('data.table')
library('mixOmics')
library('caret') # invokes xgboost
library("SuperLearner")
library('parallel')

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# loading datasets
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



#lesion <- 2-as.integer(samples$lesional)
samples$lesional_new[samples$lesional=="LES"]<- 1
samples$lesional_new[samples$lesional=="NON_LES"]<- 0   

list.keepX <- c(2:10, 15, 20)

tune <- tune.spls(X=transcriptome,
                  Y=samples$lesional_new,
                  ncomp=3,
                  test.keepX=list.keepX,
                  validation = "loo",
                  measure = "MSE",
                  progressBar = TRUE)
