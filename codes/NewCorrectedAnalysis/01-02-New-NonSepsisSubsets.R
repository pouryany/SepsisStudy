# Selecting and merging septic children expression profiles across samples.
# This code will also bring in additional samples for independent datasets. 
# The final output would be master file of children blood samples.


rm(list= ls())

library(Biobase)
library(GEOquery)
library(tidyverse)

# Loading expressions and information from original sepsis paper
# Selecting GPL570 platform only. Loading pathways table from Pathprint package

#  The following are the normal sampeles corresponding to the datasets used in 
#  Joachim et al. 2018 sepsis publication 
norm.chld           <- readRDS("DATA/normalChildrenMetaData")
colnames(norm.chld) <- c("arrayID" ,  "age", "study")
norm.chld$state     <- "normal"



# The exp.sets object contains 776 expression samples pertaining to 
# the studies used in Joachim et al (Only GPL570 platforms)
exp.sets   <- readRDS("DATA/CorrectedAllSepsisDatasets.rds")

# The final list of particular samples used in Joachim et al.
samp.tab   <- readxl::read_xlsx("DATA/TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")

# The table contains the list of Pathpring pathways and genes belonging to them.
PathExpTab <- readRDS("DATA/pathGeneTable.rds")




child.tab  <- samp.tab[samp.tab$stage == "child",]
sep.sets   <- exp.sets[,colnames(exp.sets) %in% child.tab$arrayID]
norm.exps  <- exp.sets[,colnames(exp.sets) %in% norm.chld$arrayID]



# Bringing in an additional dataset to expand control sample size

new.sets            <- getGEO("GSE35571",GSEMatrix =TRUE)
name1               <- gsub("\\_.*","",names(new.sets))
names(new.sets)     <- gsub("^.*\\.","",name1)


news.tab            <- pData(new.sets[[1]])
type.ok             <- which(news.tab$`disease status:ch1` == "non-asthmatic")
ages.ok             <- which(as.numeric(news.tab$`age (yr):ch1`)  < 11.5)
alls.ok             <- intersect(type.ok,ages.ok)
alls.ok             <- rownames(news.tab[alls.ok,])
normals.exp3        <- exprs(new.sets[[1]])[,alls.ok]


# Resolving discrepencies between gene annotations.
uni.names     <- intersect(rownames(sep.sets),rownames(normals.exp3))
sep.sets      <- sep.sets[uni.names,]
norm.exps     <- norm.exps[uni.names,]
normals.exp3  <- normals.exp3[uni.names,]


run1.exp.mat  <- do.call(cbind,list(sep.sets,norm.exps,normals.exp3))


# Processing gene ranking for pathway summary statistics.


GPL570Annotation <- readRDS("DATA/AnnotationGPL570.rds")
sym1.table       <- GPL570Annotation[,c("ID","Gene Symbol")]



sym1.table        <-  sym1.table[,c("ID","Gene Symbol")]
sym1.table        <-  sym1.table[uni.names,]
run1.exp.mat      <- as_tibble(run1.exp.mat)
run1.exp.mat$Gene <- sym1.table$`Gene Symbol`
run1.exp.mat      <- run1.exp.mat %>% group_by(.,Gene) %>%
                     summarise_if(.,is.numeric,mean)
run1.exp.mat      <- run1.exp.mat %>% filter(.,Gene != "")

run1.exp.mat      <- run1.exp.mat  %>% mutate_if(.,is.numeric,rank) 


tab.order         <- colnames(run1.exp.mat)[-c(1)]
description       <- samp.tab %>% filter(., arrayID %in% tab.order,
                                         stage == "child") %>% 
                     select(.,arrayID,age,study)
description$state <- "disease"
description       <- rbind(description,as_tibble(norm.chld))
tab.order         <- left_join(as_tibble(tab.order),
                         description, by= c("value"= "arrayID"))


saveRDS(tab.order,"DATA/CorrectedAllChildrenMeta.rds")
saveRDS(run1.exp.mat,"DATA/CorrectedAllChildren.rds")





#
#
#
#
#
#
#
#
#
#
#
#
#%>%mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

PathExpTab   <- left_join(PathExpTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))
PathExpTab   <- PathExpTab %>% group_by(., Pathway) %>% 
    summarise_if(is.numeric,mean, na.rm=T)



GPL570.tab  <- tab.order %>% group_by(.,study,state) %>% 
    summarise(., n = n())

write_csv(GPL570.tab,"Reports/NewChildGPL570samples.csv")

PathExpMat  <- as.matrix(PathExpTab[,-c(1)])

rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat <- as.data.frame(PathExpMat)

plotMA(PathExpMat)
plotDensities(PathExpMat,legend = F)
## The following is differential pathway analysis

design <- model.matrix(~ tab.order$state + 0,PathExpMat)
colnames(design) <- c("Sepsis","Normal")
fit <- lmFit(PathExpMat, design)
cont.matrix <- makeContrasts(Sepsis-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)


tT.select <- (head(tT,20))



###### Just redoing things with a Z-normalization

PathExpMat <- apply(PathExpMat, 2, function(X){ (X - mean(X))/sd(X)})
PathExpMat <- as.data.frame(PathExpMat)


#plotMA(PathExpMat2[,3:5])

design <- model.matrix(~ tab.order$state + 0,PathExpMat)
colnames(design) <- c("Sepsis","Normal")
fit <- lmFit(PathExpMat, design)
cont.matrix <- makeContrasts(Sepsis-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

head(tT,20)
hist(tT$adj.P.Val)
tT.select <- tT[ tT$adj.P.Val < 0.01 & abs(tT$logFC) > 0.4,]
nrow(tT.select)
rownames(tT.select)

write.csv(tT.select,"Reports/NewReports/DEpathwaysChildren.csv")

labs <- as.data.frame(tab.order[,2:4])
rownames(labs) <- tab.order$value
pheatmap((PathExpMat[rownames(tT.select),]),show_rownames = F,
         annotation_col = labs,clustering_method = "ward.D2")

?pheatmap

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


