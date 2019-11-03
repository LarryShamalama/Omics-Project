library(rstudioapi)
#library(gplots)
library(data.table)
library(ComplexHeatmap)
library(dendextend)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
transcriptome <- read.csv('../filter-data/transcriptome.csv')
gene_ids   <- transcriptome$X
sample_ids <- colnames(transcriptome)[-c(1)]
transcriptome <- transcriptome[,-c(1)]
transcriptome <- transpose(transcriptome)

rownames(transcriptome) <- substr(sample_ids, 1, nchar(sample_ids[1])-3)
colnames(transcriptome) <- gene_ids


genes <- read.table('../genelist.txt')$V1
colMeans(transcriptome)

transcriptome_filt <- subset(transcriptome, rownames(transcriptome) %in% genes)

euclidean <- function(x){
  return (sqrt(sum(x^2)))
}

row_dend = as.dendrogram(hclust(dist(transcriptome_filt,
                                     method="euclidean")))

Heatmap(as.matrix(transcriptome_filt),
        name='Normalized Gene\nExpression Levels',
        cluster_rows = row_dend)
