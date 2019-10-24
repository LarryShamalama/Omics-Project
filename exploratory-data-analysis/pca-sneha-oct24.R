library(rstudioapi)
library(FactoMineR)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
all_samples  <- read.csv('../filter-data/all.csv')
all_samples  <- all_samples[all_samples$clinical_group == 'AD',]

SamplesAD <- data.frame(all_samples$sample_id,
                        all_samples$clinical_group,
                        all_samples$lesional)

transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes time to load
rownames(transcriptome) <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]


transcriptome <- transcriptome[rownames(transcriptome) %in% all_samples$sample_id,]
indices <- match(all_samples$sample_id, rownames(transcriptome))
stopifnot(sum(is.na(indices)) == 0)

transcriptome <- transcriptome[indices,]
stopifnot(all(rownames(transcriptome) == all_samples$sample_id))
transcriptome <- t(transcriptome)

PCA_Transcriptome <- data.frame(t(transcriptome), lesion = SamplesAD$all_samples.lesional)
PCA_Transcriptome$lesion <- ifelse(PCA_Transcriptome$lesion == 'LES', 'Lesion', 'Non-lesion')
TransNormPCA <- PCA(PCA_Transcriptome, quali.sup = which(names(PCA_Transcriptome) == "lesion"), graph=FALSE)

plot(TransNormPCA, 
     choix = "ind", 
     label = "none", 
     habillage = which(names(PCA_Transcriptome) == "lesion"))
