library(rstudioapi)
library(pheatmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
transcriptome <- read.csv('../filter-data/transcriptome.csv')
rownames(transcriptome) <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]

#pheatmap(transcriptome) #takes too much time