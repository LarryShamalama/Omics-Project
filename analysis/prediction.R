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

n_components <- 10
n_keep       <- 10

spls <- mixOmics::splsda(transcriptome, 
                         samples$lesional_new, 
                         ncomp=n_components,
                         keepX=rep(n_keep, n_components))

cvspls <- perf(spls, validation="loo")

cumul_gene_names <- c()

for (i in 1:n_components){
  temp_component <- spls$loadings$X[,i]
  temp_component <- temp_component[temp_component != 0]
  cumul_gene_names <- c(names(temp_component), cumul_gene_names)
}

# list of genes selected by spls
cumul_gene_names <- unique(cumul_gene_names) 

write.table(cumul_gene_names, file='../genelist_spls.txt')

cat('Number of selected genes:', length(cumul_gene_names), '\n')
cat('Maximum number of selected genes:', n_components*n_keep, '\n')

models <- list('SL.ranger', 
               'SL.glm', 
               #'SL.ksvm',
               'SL.xgboost',
               'SL.nnet',
               'SL.mean')


sl <- SuperLearner(Y=lesion.status,
                   X=transcriptome[cumul_gene_names],
                   family=binomial(),
                   SL.library=models)

if (parallel::detectCores() > 8){
  num_folds <- length(lesion.status)
  cat('Using multiprocessing for leave-one-out cross-validation\n')
} else{
  num_folds <- 5
  cat('Not enough computational power... using 5-fold cross-validation')
}

cv.sl <- CV.SuperLearner(Y=lesion.status,
                         X=transcriptome[cumul_gene_names],
                         V=num_folds, # ideally use n for leave-out-one cross-validation
                         family=binomial(),
                         SL.library=models)

predictions <- predict.SuperLearner(cv.sl)$pred

source('utils.R') # importing functions from another R script

plot.roc(predictions, lesion.status)
ggsave('roc_curve.png', width=7, height=4.5, units='in', dpi=250)