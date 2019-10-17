

# set directory
setwd("/Users/thewooz/Documents/Omics-Project/data")

# import data
gene_samples<- read.csv("M2PHDS_19-20_OMICS_TRANSC_MAARS_normTranscriptome_618samples_16042014.txt", sep="\t")

# create a random sample of the samples
randomsamples<- gene_samples[ ,sample(ncol(gene_samples), 10) ]
# summarize the random samples
summary(randomsamples)


#variance
V <- apply(randomsamples, 1, var)

#order by decreasing variance and place the names of those genes in a vector
selectedGenes <- names(V[order(V, decreasing = T)][1:100])


library(pheatmap)
pheatmap(randomsamples[selectedGenes,], scale = 'row', show_rownames = FALSE)


# boxplot
boxplot(gene_samples)

, xlab="Samples", ylab="log2 (raw counts)",axes=FALSE)
axis(2)
box()
legend("topright", c("London","South Africa","Leicester"), col=c("blue2","gold","green2"), pch=15, horiz=TRUE, bg=NA)


## DENDOGRAM
deucl_raw<- dist(t(select(expr_singhania_raw_no0_cpm2p_prot, -starts_with("Gene"))))
deucl_norm<- dist(t(expr_singhania_raw_no0_cpm2p_prot_TMMlog2cpm))

plot(flashClust::hclust(deucl_raw, method="ward"), labels=infos_singhania$Location)
plot(flashClust::hclust(deucl_norm, method="ward"), labels=infos_singhania$Location)

