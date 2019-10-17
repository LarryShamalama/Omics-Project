# Omics Presentation - Data Exploration 
#  K. Creel
#  17 Oct 2019


library(ggplot2)

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read in the cleaned data
transcriptome<- read.csv("../filter-data/transcriptome.csv")
description<- read.csv("../filter-data/all.csv")

# check for NAs in the transcriptome data
NAs<- sum(is.na(transcriptome))
NAs

# order the clinical groups
description$clinical_group<- factor(description$clinical_group, levels=c("AD","PSO","CTRL"))


##  DESCRIBE THE SAMPLES

# table of percent of SAMPLES per clinical group and lesional/non-lesional
prop.table(table(description$clinical_group,description$lesional))*100

# bar chart for SAMPLES per clinical group and lesional/non-lesionsl
ggplot2::ggplot(data=description, aes(clinical_group, fill=lesional))+
  geom_bar(stat="count", position="stack")+
  scale_fill_manual(values=c("red4","steelblue4"),name="Lesioned/Non-Lesioned", labels=c("Lesioned","Non-Lesioned"))+
  theme_minimal()+
  xlab("Clinical Group")+
  ylab("Number of Samples")

 
##  VERIFY NORMALITY

bp_col<- description$lesional
levels(bp_col) <- c("red4", "steelblue4")
boxplot(transcriptome[,-1], col=as.character(bp_col), 
        xlab="Samples", ylab="Nomalized CPM", axes=F)
axis(2)
box()
#legend("topright", c("Lesioned", "Non-Lesioned"), col = c("red4","steelblue4"), pch = 15, horiz = TRUE, bg = NA)

# dendogram
dendo_data <- dist(t(transcriptome[,-1]))
plot(flashClust::hclust(dendo_data, method = "ward"), labels = description$lesioned)







