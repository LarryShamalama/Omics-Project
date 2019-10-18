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
    BiocManager::install(pkg)
  })
}

tryCatch({
  prin
}, error=function(e){
  cat('hello')
})
library('ggplot2')
library('rstudioapi')
library('data.table')
library('mixOmics')
library('caret') # invokes xgboost
library("SuperLearner")
library('parallel')

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# loading datasets
all_patients  <- read.csv('../filter-data/all.csv')
all_samples   <- read.csv('../filter-data/ad_full.csv')
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

cat('Number of selected genes:', length(cumul_gene_names), '\n')
cat('Maximum number of selected genes:', n_components*n_keep, '\n')

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
                         SL.library=models,
                         parallel='multicore')

predictions <- predict.SuperLearner(cv.sl)$pred

source('utils.R') # importing functions from another R script

plot.roc(predictions, lesion.status)
ggsave('roc_curve.png', width=7, height=4.5, units='in', dpi=250)