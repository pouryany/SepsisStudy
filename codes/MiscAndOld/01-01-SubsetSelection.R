# Just a separation of data based on the platform. 
# The current study mostly focuses on GPL570 Platform
# By the end of this code, I will have a table gene expressions.

rm(list= ls())

library(Biobase)
library(tidyverse)
library(limma)

exp.sets <- readRDS("DATA/AllSepsisDatasets.rds")
samp.tab <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")


exp.sets[1]

exp.normalizer <- function(gset){
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }
    return(gset)
}

exp.sets       <- lapply(exp.sets, exp.normalizer)

names(exp.sets)

#fData.run1 <- fData(exp.sets[[1]])
#fData.run1$ID
exp.mat.sets <- sapply(exp.sets,function(X){ (exprs(X))})
#head(exp.mat.sets)
#sapply(exp.mat.sets,nrow)
name1 <- gsub("\\_.*","",names(exp.mat.sets))
names(exp.mat.sets) <- gsub("^.*\\.","",name1)
names(exp.sets)     <- gsub("^.*\\.","",name1)
#samp.tab %>% filter(., platform == "GPL571") %>% select(.,study) %>% unique()
# The above list contains all expression matrices from the study

# Lets focus only on the GPL570 Samples, hereby, run1


run1.list    <-  samp.tab %>% filter(., platform == "GPL570")
fData.run1   <-  fData(exp.sets[[unique(run1.list$study)[1]]])

saveRDS(fData.run1,"DATA/AnnotationGPL570.rds")

sym1.table   <-  fData.run1[,c("ID","Gene Symbol")]
run1.exp.mat <-  sapply(exp.mat.sets, function(X){
                        X[,colnames(X) %in% run1.list$arrayID]})
run1.exp.mat <- run1.exp.mat[sapply(run1.exp.mat,ncol) != 0]


# Testing for normalizations 
    names(run1.exp.mat)
    test.set <- run1.exp.mat[[1]]
    ex <- (test.set)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    test.set <- log2(ex) }
    
    limma::plotMA((test.set))
    limma::plotDensities((test.set),legend = F)
    
 
## All data from Wong needs reinvestigation, maybe?  
       
    
run1.exp.mat <- do.call(cbind, run1.exp.mat)


run2.list    <-  samp.tab %>% filter(., platform == "GPL571")
fData.run2   <-  fData(exp.sets[[unique(run2.list$study)[1]]])
sym2.table   <-  fData.run2[,c("ID","Gene Symbol")]
run2.exp.mat <-  sapply(exp.mat.sets, function(X){
                       X[,colnames(X) %in% run2.list$arrayID]})
run2.exp.mat <- run2.exp.mat[sapply(run2.exp.mat,ncol) != 0]
run2.exp.mat <- do.call(cbind, run2.exp.mat)

#rm(exp.sets)

#sum(rownames(run1.exp.mat) != sym1.table$ID)
#sum(rownames(run2.exp.mat) != sym2.table$ID)

run1.exp.mat      <- as_tibble(run1.exp.mat)
run1.exp.mat$Gene <- sym1.table$`Gene Symbol`
run1.exp.mat      <- run1.exp.mat %>% group_by(.,Gene) %>%
                     summarise_if(.,is.numeric,mean)
run1.exp.mat      <- run1.exp.mat %>% filter(.,Gene != "")


run2.exp.mat      <- as_tibble(run2.exp.mat)
run2.exp.mat$Gene <- sym2.table$`Gene Symbol`
run2.exp.mat      <- run2.exp.mat %>% group_by(.,Gene) %>%
                     summarise_if(.,is.numeric,mean)
run2.exp.mat      <- run2.exp.mat %>% filter(.,Gene != "")


saveRDS(run1.exp.mat,"DATA/SepsisGPL570.rds")
saveRDS(run2.exp.mat,"DATA/SepsisGPL571.rds")


##
#
#
#
#run1.exp.mat <- readRDS("run1.rds")