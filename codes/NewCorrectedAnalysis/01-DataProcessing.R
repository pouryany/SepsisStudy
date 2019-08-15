rm(list = ls())
gc()


library(Biobase)
library(tidyverse)
library(limma)


# Just a separation of data based on the platform. 
# The current study mostly focuses on GPL570 Platform
# By the end of this code, I will have a table gene expressions.


# The below datasets are precalculated. 
# TO BE UPDATED.  REFERENCING THE CODES TO GET THESE THINGS. 
# DOCUMENTING THE DATA

exp.mat          <- readRDS("DATA/CorrectedAllSepsisDatasets.rds")
GPL570Annotation <- readRDS("DATA/AnnotationGPL570.rds")
sym1.table       <- GPL570Annotation[,c("ID","Gene Symbol")]
samp.tab         <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")



# Getting only the sepsis expression profiles
run1.list        <-  samp.tab %>% filter(., platform == "GPL570")
sym1.table       <-  sym1.table[,c("ID","Gene Symbol")]
run1.exp.mat     <-  exp.mat[,colnames(exp.mat) %in% run1.list$arrayID]



run1.exp.mat      <- as_tibble(run1.exp.mat)
run1.exp.mat$Gene <- sym1.table$`Gene Symbol`
run1.exp.mat      <- run1.exp.mat %>% group_by(.,Gene) %>%
                     summarise_if(.,is.numeric,mean)
run1.exp.mat      <- run1.exp.mat %>% filter(.,Gene != "")

saveRDS(run1.exp.mat,"DATA/CorrectedSepsisGPL570.rds")

