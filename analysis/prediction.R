# reset console
cat('\014')
rm(list=ls())

library('ggplot2')
library('rstudioapi')
library('data.table')
library('mixOmics')
library('caret') # invokes xgboost
library("SuperLearner")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# loading datasets
all_patients  <- read.csv('../filter-data/all_patients.csv')
all_samples   <- read.csv('../filter-data/all_samples.csv')
transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes a while to load

stopifnot(dim(transcriptome)[1] == dim(all_samples)[1])

sample_id <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]
rownames(transcriptome) <- sample_id

lesion.status <- all_samples$lesional
lesion.status <- as.integer(ifelse(lesion.status == 'LES', 1, 0)) # binary encoding

n_components <- 8
n_keep       <- 4

spls <- mixOmics::splsda(transcriptome, 
                         lesion.status, 
                         ncomp=n_components,
                         keepX=rep(n_keep, n_components))

cumul_gene_names <- c()

for (i in 1:n_components){
  temp_component <- spls$loadings$X[,i]
  temp_component <- temp_component[temp_component != 0]
  cumul_gene_names <- c(names(temp_component), cumul_gene_names)
}

cumul_gene_names <- unique(cumul_gene_names)

cat('Number of selected genes:', length(cumul_gene_names))
cat('Maximum number of selected genes:', n_components*n_keep)

models <- list('SL.ranger', 
               'SL.glm', 
               'SL.ksvm',
               'SL.xgboost',
               'SL.nnet',
               'SL.mean')

sl <- SuperLearner(Y=lesion.status,
                   X=transcriptome[cumul_gene_names],
                   family=binomial(),
                   SL.library=models)

cv.sl <- CV.SuperLearner(Y=lesion.status,
                         X=transcriptome[cumul_gene_names],
                         V=5,
                         family=binomial(),
                         SL.library=models)
