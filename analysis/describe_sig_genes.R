

# select the transcriptomes for the significant genes
gene_names_at<- rownames(gene_pvalues)
trans_sig <- transcriptome[rownames(transcriptome) %in% gene_names_at,]

# dendrogram of significant genes
dendo_data <- dist(t(trans_sig[,-1]))
plot(flashClust::hclust(dendo_data, method = "ward"))


dotchart(exprs(rownames(transcriptome)))
         