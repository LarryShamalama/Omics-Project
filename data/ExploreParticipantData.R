
library("ggplot2")

# set directory
setwd("/Users/thewooz/Documents/Omics-Project/data")

# import data
features_all_samples<- read.csv("M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_all_Fri_Apr_04_14h_CEST_2014.csv", sep="\t")
features_AD_samples<- read.csv("M2PHDS_19-20_OMICS_CLIN_DATA_MAARS_AD_full_20190131_12-34-49.csv", sep="\t")
gene_samples<- read.csv("M2PHDS_19-20_OMICS_TRANSC_MAARS_normTranscriptome_618samples_16042014.txt", sep="\t")



## EXPLORE PARTICIPANTS

# create dataset of unique participants
reduce<- features_all_samples[ ,c(6,2,7:34)]
undup_pts<- reduce[!duplicated(reduce$MAARS_identifier), ]

# AD and control group
addmargins(table(undup_pts$clinical_group))

# gender (percent)
round(prop.table(table(undup_pts$Gender))*100, 2)

# age (mean)
round(mean(undup_pts$CUSTOM_Age, na.rm=TRUE))

# age (summary/spread)
summary(undup_pts$CUSTOM_Age)

# institution (percent)
round(prop.table(table(undup_pts$Institution))*100, 2)

# Global Assessment Score
table(undup_pts$Global_Assessment_Score)
round(prop.table(table(undup_pts$Global_Assessment_Score))*100,2)

#  plot of GAS per clinical group
ggplot(undup_pts) +
  geom_bar(aes(x=clinical_group, fill=Global_Assessment_Score), position="fill") +
  scale_x_discrete(limits=c("AD","PSO","CTRL")) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_brewer()+
  theme(panel.background = element_rect(fill = "grey95"), 
        panel.grid.major = element_line(colour = "grey70"), 
        plot.background  = element_rect(fill="white")) +
  ylab("Clinical Group") +
  xlab("Percent")


# Known allergies
library(qdapTools)
allergies<- (mtabulate(undup_pts[grep('Allergies',names(undup_pts), value=T)]))
allergies$percent<- allergies[,2]/allergies[,1]+allergies[,2]


# Medications
meds<- (mtabulate(undup_pts[grep('Medication',names(undup_pts), value=T)]))
meds$percent<- meds[,2]/meds[,1]+meds[,2]


# Chronic Disease
chrondis<- (mtabulate(undup_pts[grep('chronic',names(undup_pts), value=T)]))
chrondis$percent<- chrondis[,2]/chrondis[,1]+chrondis[,2]


# Chronic Disease
malig<- (mtabulate(undup_pts[grep('Malig',names(undup_pts), value=T)]))
malig$percent<- malig[,2]/malig[,1]+malig[,2]


# Family History of AD or Psoriasis
FamHist<- (mtabulate(undup_pts[grep('Fam',names(undup_pts), value=T)]))
FamHist$percent<- FamHist[,2]/FamHist[,1]+FamHist[,2]


### How to group the above tables by clinical_group? Getting errors when piping...
library("dplyr")
FamHist<- undup_pts %>%
          group_by(undup_pts$clinical_group) %>%
          mtabulate(undup_pts[grep('Fam',names(undup_pts), value=T)])

