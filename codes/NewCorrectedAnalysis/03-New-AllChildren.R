# Just a separation of data based on the platform. 
# The current study mostly focuses on GPL570 Platform

rm(list= ls())

library(Biobase)
library(GEOquery)
library(tidyverse)
library(pheatmap)

tab.order    <- readRDS("DATA/CorrectedAllChildrenMeta.rds")
samp.tab     <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")
left.out     <- samp.tab %>% filter(.,study == "GSE26378") %>% select(.,arrayID) %>% pull()

run1.exp.mat <- readRDS("DATA/CorrectedAllChildren.rds")
PathExpTab   <- readRDS("DATA/pathGeneTable.rds")


run1.exp.mat <- run1.exp.mat  %>%  mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

PathExpTab   <- left_join(PathExpTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))
PathExpTab   <- PathExpTab %>% group_by(., Pathway) %>% 
                summarise_if(is.numeric,mean, na.rm=T)


PathExpMat  <- as.matrix(PathExpTab[,-c(1)])

rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat <- as.data.frame(PathExpMat)

plotMA(PathExpMat)
plotDensities(PathExpMat,legend = F)
## The following is differential pathway analysis

# design <- model.matrix(~ tab.order$state + 0,PathExpMat)
# colnames(design) <- c("Sepsis","Normal")
# fit <- lmFit(PathExpMat, design)
# cont.matrix <- makeContrasts(Sepsis-Normal, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# 

#tT.select <- (head(tT,50))



###### Just redoing things with a Z-normalization

PathExpMat <- apply(PathExpMat, 2, function(X){ (X - mean(X))/sd(X)})
PathExpMat <- as.data.frame(PathExpMat)


left.out1   <- setdiff(colnames(PathExpMat), left.out)
PathExpMat <- PathExpMat[,left.out1]

tab.order  <-tab.order[tab.order$value %in% left.out1,]


#plotMA(PathExpMat2[,3:5])

design <- model.matrix(~ tab.order$state + 0,PathExpMat)
colnames(design) <- c("Sepsis","Normal")
fit <- lmFit(PathExpMat, design)
cont.matrix <- makeContrasts(Sepsis-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT   <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

head(tT,20)
hist(tT$adj.P.Val)

{
  adt.chd <- read.csv("Reports/NewReports/DEpathways.csv")
  adt.chd <- tT[as.character(adt.chd$X)[4],]  
  nrow(adt.chd[adt.chd$adj.P.Val < 0.05 & abs(adt.chd$logFC) > 0.1,])
}


tT.select <- tT[ tT$adj.P.Val < 0.001 & abs(tT$logFC) > 0.3,]
nrow(tT.select)
rownames(tT.select)

tab.order <- mutate(tab.order, age = as.numeric(age))
test <- left_join(tab.order,samp.tab[,c(1,6)],by = c("value" ="arrayID"))
labs <- as.data.frame(test[,2:4])

rownames(labs) <- tab.order$value
pheatmap((PathExpMat[rownames(tT.select),]),show_rownames = F,show_colnames = F,
         annotation_col = labs,clustering_method = "ward.D2",scale = "row")




#Subtypes exist in children. Not much in adults. 


datasets <- samp.tab %>% filter(., platform == "GPL570") %>% 
                          select(.,study) %>% unique() %>% pull

name1    <- gsub("\\_.*","",names(exp.sets))
names(exp.sets)     <- gsub("^.*\\.","",name1)

exp.sets <- exp.sets[datasets]

# Getting an additional dataset to boost the number of normals



# GSE26378 contains age. Disease state which contains normal control
# GSE26440 contains age. Disease state which contains normal control
norm.sets <- exp.sets[1:2]



# Looking at online datasets for whole blood samples
# GSE35571: Childern Asthma/Non-asthma with age variable. Many samples
# GSE58667: Control Samples with age
# GSE14844: Many 12 year olds here
# GSE72439: A lot of samples of patients with Growth hormone deficiency
# GSE6575: Normal samples of children but no age variable on first glance
# GSE11755: Contains children sepsis, time series, and normal, no age variable
# GSE8121: Hector Wong dataset
# GSE67684 not related but look at this dataset...
#
# GSE9692, GSE13904, and GSE4607 have median age below 5. The age data exists
# but is not available online. Maybe Contact for them is Hector Wong?


