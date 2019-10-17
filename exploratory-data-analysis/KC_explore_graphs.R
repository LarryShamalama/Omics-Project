
library("ggplot2")
library('rstudioapi')
library('dplyr')
library('flashClust')

# set directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# import data
features_all_samples<- read.csv("M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_all_Fri_Apr_04_14h_CEST_2014.csv", sep="\t")
features_AD_samples<- read.csv("M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_AD_full_20190131_12-34-49.csv", sep="\t")
gene_samples<- read.csv("M2PHDS_19-20_OMICS_TRANSC_MAARS_normTranscriptome_618samples_16042014.txt", sep="\t")



# verify normality
#  COLORS PER 'lesional' ARE NOT WORKING

bp_col<- as.factor(=features_AD_samples$lesional)
levels(bp_col) <- c("blue2", "gold", "green2")
boxplot(gene_samples, col=as.character(bp_col), 
        xlab="Samples", ylab="Nomalized CPM", axes=F)
axis(2)
box()



# dendogram
#  COLORS NOT WORKING
dendo<- dist(t(gene_samples))
plot(flashClust::hclust(dendo, method = "ward"), labels = features_all_samples$lesional)

#Boris' code
#  deucl_norm <- dist(t(expr_singhania_no0_cpm2p_prot_TMMlog2cpm))
#  plot(flashClust::hclust(deucl_norm, method = "ward"), labels = infos_singhania$Location)





