
rm(list = ls())

library(pheatmap)
library(limma)
library(tidyverse)

run1.exp.mat   <- readRDS("DATA/normalChildrenExpression.rds")
run1.exp.mat   <- run1.exp.mat  %>% mutate_if(.,is.numeric,rank) 
labs           <- readRDS("data/normalChildrenMetaData")
rownames(labs) <- labs$ArrayID
labs$ArrayID   <- NULL
labs$Age       <- as.numeric(as.character(labs$Age))
#%>%mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

PathExpTab  <- readRDS("DATA/pathGeneTable.rds")
PathExpTab  <- left_join(PathExpTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))
PathExpTab  <- PathExpTab %>% group_by(., Pathway) %>% 
               summarise_if(is.numeric,mean, na.rm=T)


PathExpMat  <- as.matrix(PathExpTab[,2:ncol(PathExpTab)])

rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat2 <- as.data.frame(PathExpMat)
PathExpMat2 <- apply(PathExpMat2, 2, function(X){ (X - mean(X))/sd(X)})
PathExpMat2 <- as.data.frame(PathExpMat2)



selectpaths <- read.csv("Reports/S3_LassoChildrenPathways.csv")

pheatmap(PathExpMat2[as.character(selectpaths$Pathway.Name),],show_rownames = F, clustering_method = "ward.D2",
         annotation_col = labs )


plotDensities(apply(PathExpMat2, 1, var))
