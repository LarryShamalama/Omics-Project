# reset console
cat('\014')
rm(list=ls())

library('ggplot2')
library('rstudioapi')
library('data.table')
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# loading datasets
all_patients <- read.csv('../filter-data/all_patients.csv')
all_samples  <- read.csv('../filter-data/all_samples.csv')
transcriptome <- read.csv('../filter-data/transcriptome.csv') # takes a while to load

