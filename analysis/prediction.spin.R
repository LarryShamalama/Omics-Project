# reset console
cat('\014')
rm(list=ls())

set.seed(1)
selected_genes <- c('ENSG00000186832', 
                    'ENSG00000124102', 
                    'ENSG00000165474', 
                    'ENSG00000172382', 
                    'ENSG00000198074', 
                    'ENSG00000163220', 
                    'ENSG00000163221') # manually entered

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
colnames(transcriptome) <- substr(colnames(transcriptome), 1, nchar(colnames(transcriptome)[1])-3)
transcriptome <- transcriptome[,colnames(transcriptome) %in% selected_genes]
# stop if not every row matches a sample_id
stopifnot(all(rownames(transcriptome) == samples$sample_id))



#lesion <- 2-as.integer(samples$lesional)
samples$lesional_new[samples$lesional=="LES"]<- 1
samples$lesional_new[samples$lesional=="NON_LES"]<- 0   


models <- list('SL.ranger', 
               'SL.glm', 
               'SL.ksvm',
               'SL.xgboost',
               'SL.nnet',
               'SL.mean')

if (parallel::detectCores() > 8){
  num_folds <- length(lesion.status)
  cat('Using multiprocessing for leave-one-out cross-validation\n')
} else{
  num_folds <- 10
  cat('Not enough computational power... using 5-fold cross-validation')
}


sl <- SuperLearner(Y=samples$lesional_new,
                   X=transcriptome,
                   family=binomial(),
                   SL.library=models)

loo.sl <- CV.SuperLearner(Y=samples$lesional_new,
                          X=transcriptome,
                          family=binomial(),
                          V=length(samples$lesional_new),
                          verbose=TRUE,
                          SL.library=models)

bootstrap <- function(X, Y){
  resample <- sample(1:length(Y), replace=TRUE)
  
  return (list(X[resample,], Y[resample]))
}

predictions  <- predict.SuperLearner(sl)$pred
num_resample <- 10

bs_predictions <- c()

for (i in 1:num_resample){
  bs_sample <- bootstrap(transcriptome, samples$lesional_new)
  X_temp <- bs_sample[[1]]
  Y_temp <- bs_sample[[2]]
  temp.sl <- SuperLearner(Y=Y_temp,
                          X=X_temp,
                          family=binomial(),
                          SL.library=models)
  
  bs_predictions <- cbind(predict.SuperLearner(temp.sl)$pred,
                          bs_predictions)
  if (i %% 10 == 0){
    cat('Done', i, 'iterations\n')
  }
}

source('utils.R') # importing functions from another R script

plot.roc(predict.SuperLearner(sl)$pred, 
                  samples$lesional_new)
ggsave('roc_curve.png', width=7, height=4.5, units='in', dpi=250)