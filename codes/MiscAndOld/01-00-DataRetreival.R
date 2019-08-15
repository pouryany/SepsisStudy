rm(list = ls())

library(Biobase)
library(GEOquery)
library(limma)
library(tidyverse)

samp.tab <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")
datasets <- unique(samp.tab$study)
# The next three lines are for generating data. Takes time
exp.sets <- sapply(datasets,getGEO,GSEMatrix =TRUE)
exp.sets <- unlist(exp.sets)
saveRDS(exp.sets,"DATA/AllSepsisDatasets.rds")
