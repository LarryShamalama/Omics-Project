# reset console
cat('\014')

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
list.keepX <- c(2:10, 15, 20)

set.seed(1)

tune.Mfold <- tune.spls(X=transcriptome,
                        Y=samples$lesional_new,
                        ncomp=n_components,
                        test.keepX=list.keepX,
                        validation = "Mfold",
                        folds=10, # initially 10
                        measure = "MSE",
                        progressBar = TRUE)


cumul_gene_names <- c()

spls <- mixOmics::splsda(transcriptome, 
                         samples$lesional, 
                         ncomp=n_components,
                         keepX=as.vector(tune.Mfold$choice.keepX))

plotIndiv(spls, comp=1:2, ind.names=FALSE, rep.space='X-variate', legend=TRUE)
plotVar(spls)

for (i in 1:n_components){
  temp_component <- spls$loadings$X[,i]
  temp_component <- temp_component[temp_component != 0]
  cumul_gene_names <- c(names(temp_component), cumul_gene_names)
}

write.table(cumul_gene_names, file='../genelist_spls_Mfold.txt')
