#  Gene Differentiation 
#  L.Dong and K.Creel
#  24-10-2019

# clear the environment
# clear the console
cat('\014')


# load libraries
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(rstudioapi)
library(limma)
library(data.table)

# read in sample information, keeping those that are AD
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
# transpose transcriptome data
transcriptome <- t(transcriptome)


#lesion <- 2-as.integer(samples$lesional)
samples$lesional_new[samples$lesional=="LES"]<- 1
samples$lesional_new[samples$lesional=="NON_LES"]<- 0   

library(gdata)
# create the model matrix for paired data
#   drop the unused dummy variables with drop.levels
paired.design<- as.data.frame(model.matrix(~lesional_new + drop.levels(MAARS_identifier), data=samples))

# fitting linear model
fit <- lmFit(transcriptome, design=paired.design)
# squeezes genewise residual variance towards a common value
fit <- eBayes(fit, robust=TRUE)

# significance tests - identify which genes are differentially expressed
signif <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(signif)[ ,2]


# extract genes p-values 
#  NOTE:  code is giving ONLY the significant genes
top.table <- topTable(fit, adjust.method="BH", p.value = 0.05, lfc = 1, sort.by = "P", n=Inf, coef=2)
top.table2 <- topTable(fit, adjust.method="BH", p.value = 1, lfc = 0, sort.by = "P", n=Inf, coef=2)


# plotting histograms
library(ggplot2)
n <- length(top.table2$adj.P.Val)
df <- data.frame(cbind(c(top.table2$P.Value, top.table2$adj.P.Val), (c(rep('p-value', n), rep('adj. p-value', n)))))
colnames(df) <- c('p_value', 'cat')

ggplot(df, aes(x=cat, y=p_value)) +
  geom_boxplot()


gene_pvalues_all<- top.table[ ,(ncol(top.table)-3):ncol(top.table)]
#gene_pvalues_sig<- gene_pvalues_all[which(gene_pvalues_all$adj.P.Val<.05), ]
gene_pvalues_sig<- top.table[ ,(ncol(top.table)-3):ncol(top.table)]

# Arul's code:
topTable(BayesFit, adjust.method = "BH")

test<- topTable(fit, adjust.method = "BH", n=Inf)
test2<- test[ ,(ncol(test)-3):ncol(test)]
names(test)
summary(test)

# check to see how many genes are p<.05
length(which(test2$P.Val < 0.05))


# save gene pvalues
save(gene_pvalues,file="../filter-data/gene_pvalues.Rda")


# remove the _at at the end of gene (only do this when exporting, otherwise the gene name won't match the other files)
gene_names<- substr(row.names(gene_pvalues),1,nchar(row.names(gene_pvalues))-3)

# save gene names to a .csv
write.csv(rownames(gene_pvalues), "/Users/thewooz/Documents/Omics-Project/filter-data/genelist.csv")

