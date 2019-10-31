

# set directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# import data
#gene_samples<- read.csv("M2PHDS_19-20_OMICS_TRANSC_MAARS_normTranscriptome_618samples_16042014.txt", sep="\t")

# read in sample information, keeping those that are AD
samples <- read.csv('../filter-data/all.csv')
samples_AD <- subset(samples, samples$clinical_group == 'AD')

# read in transcriptome data, change row names and remove the var 'X'
transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes time to load
rownames(transcriptome) <- transcriptome$X
transcriptome <- transcriptome[,-c(1)]

# keep the transcriptomes that are AD
transcriptome_AD <- transcriptome[rownames(transcriptome) %in% samples_AD$sample_id,]
# create an index, matching the sample_id to the row in the transcriptome
indices <- match(samples_AD$sample_id, rownames(transcriptome_AD))
# stop the process if there's any lines that have NA (sample_id has not matching transcriptome)
stopifnot(sum(is.na(indices)) == 0)

# order the transcriptome data according to the vector
transcriptome_AD <- transcriptome_AD[indices,]
# stop if not every row matches a sample_id
stopifnot(all(rownames(transcriptome_AD) == samples_AD$sample_id))
# transpose transcriptome data
transcriptome_AD <- t(transcriptome_AD)




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



## DENDOGRAM
dendo<- dist(transcriptome_AD)

plot(flashClust::hclust(deucl_raw, method="ward"), labels=infos_singhania$Location)
plot(flashClust::hclust(deucl_norm, method="ward"), labels=infos_singhania$Location)
