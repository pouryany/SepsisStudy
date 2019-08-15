# This is the code for downloading the data in the sepsis paper. 
# The master file for the samples was provided by Les, TrueSepsisData...
# This code downloads all of the datasets. So, watch out for unwanted samples.


rm(list = ls())
# The lists of datasets in this analaysis



            # gset <- getGEO("GSE28750", GSEMatrix =TRUE, AnnotGPL=TRUE)
            # if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
            # gset <- gset[[idx]]
            # fvarLabels(gset) <- make.names(fvarLabels(gset))
            # 
            # # log2 transform
            # ex <- exprs(gset)
            # qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
            # LogC <- (qx[5] > 100) ||
            #     (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            #     (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
            # if (LogC) { ex[which(ex <= 0)] <- NaN
            # exprs(gset) <- log2(ex) }
            # 
            # 
            # 
            # gsms1 <- "00000000001111111111111111111122222222222"
            # sml1  <- unlist(str_split(gsms1,""))
            # sml1  <- gsub("0","S",x = sml1)
            # sml1  <- gsub("1","N",x = sml1)
            # sml1  <- gsub("0","O",x = sml1)
            # 
            # f.this <- fData(gset)
            # head(exprs(gset))


library(Biobase)
#BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
library(tidyverse)

samp.tab <- readxl::read_xlsx("TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")
datasets <- unique(samp.tab$study)
# The next three lines are for generating data. Takes time
#exp.sets <- sapply(datasets,getGEO,GSEMatrix =TRUE)
#exp.sets <- unlist(exp.sets)
#saveRDS(exp.sets,"DATA/AllSepsisDatasets.rds")

exp.sets <- readRDS("DATA/AllSepsisDatasets.rds")
        #this.names <- samples.tab$arrayID[samples.tab$study == datasets[1]]
        #attributes(test.this.set)
        #test.this.set <- this.set[1]
        #xx <- attr(test.this.set,"names")
        #exprs((test.this.set[[xx]]))

fData.run1 <- fData(exp.sets[[1]])
fData.run1$
exp.mat.sets <- sapply(exp.sets,function(X){ (exprs(X))})
head(exp.mat.sets)


    sapply(exp.mat.sets,nrow)
    samp.tab %>% filter(., platform == "GPL571") %>% select(.,study) %>% unique()
# The above list contains all expression matrices from the study

# Lets focus only on the GPL570 Samples, hereby, run1

run1.list  <-  samp.tab %>% filter(., platform == "GPL570")
run1.list  <-  run1.list$arrayID 
fData.run1 <- fData(exp.sets[["GSE28750.GSE28750_series_matrix.txt.gz"]])
sym.table  <- fData.run1[,c("ID","Gene Symbol")]

rm(exp.sets)


run1.exp.mat <- sapply(exp.mat.sets, function(X){
                       X[,colnames(X) %in% run1.list]})
run1.exp.mat <- run1.exp.mat[sapply(run1.exp.mat,ncol) != 0]
run1.exp.mat <- do.call(cbind, run1.exp.mat)


ncol(run1.exp.mat)      

library(limma)        

        #which(run1.exp.mat >20,arr.ind = T)
        #((run1.exp.mat[,80]))

# There is something wrong here. I cannot tell... Maybe ask John H.
# The density plot of the probes looks dramatically different
# The cut-off 100000 on the sum of the expressions seems to separate 


seems.ok  <- colnames(run1.exp.mat[,colSums(run1.exp.mat) > 100000])
seems.ok2 <-samp.tab %>% filter(.,arrayID %in% seems.ok) 
# Two studies containing adults seem ok on probe distribution


seems.dif  <- colnames(run1.exp.mat[,colSums(run1.exp.mat) < 100000])
seems.dif2 <-samp.tab %>% filter(.,arrayID %in% seems.dif) 
# All studies in this different categories are from the childrenn





# The below line shows if the two lists are aligned

#run1.exp.mat <- readRDS("run1.rds")
sum(rownames(run1.exp.mat) != sym.table$ID)

run1.exp.mat      <- as_tibble(run1.exp.mat)
run1.exp.mat$Gene <- sym.table$`Gene Symbol`
run1.exp.mat      <- run1.exp.mat %>% group_by(.,Gene) %>%
                     summarise_if(.,is.numeric,mean)
run1.exp.mat      <- run1.exp.mat %>% filter(.,Gene != "")

#saveRDS(run1.exp.mat,"DATA/run1.rds")
# Start Here For Now
run1.exp.mat <- readRDS("DATA/run1.rds")

run1.exp.mat <- run1.exp.mat  %>% mutate_if(.,is.numeric,rank) 
                #%>%mutate_if(.,is.numeric,function(X){X*X})
                # Edit ranking scales above!

PathGeneTab  <- readRDS("DATA/pathGeneTable.rds")
pathGeneExp  <- left_join(PathGeneTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))


PathExpTab  <- pathGeneExp %>% group_by(., Pathway) %>% 
               summarise_if(is.numeric,mean, na.rm=T)


plotDensities(PathExpTab[,2:100],legend = F)

tab.order   <- colnames(PathExpTab)
tab.order   <- tab.order[-c(1)]
description <- samp.tab %>% filter(., arrayID %in% tab.order) %>% 
               select(.,arrayID,stage,study,outcome)

tab.order2  <- left_join(as_tibble(tab.order), description, by= c("value"= "arrayID"))
tab.order   == tab.order2$value

PathExpMat <- as.matrix(PathExpTab[,2:132])

rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat <- as.data.frame(PathExpMat)

design <- model.matrix(~ tab.order2$stage + 0,PathExpMat)
colnames(design) <- c("AD","CH")
fit <- lmFit(PathExpMat, design)
cont.matrix <- makeContrasts(AD-CH, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)


head(tT)



###### Just redoing things with a Z-normalization

    PathExpMat2 <- as.data.frame(PathExpMat)
    PathExpMat2 <- apply(PathExpMat2, 2, function(X){ (X - mean(X))/sd(X)})
    PathExpMat2 <- as.data.frame(PathExpMat2)
    design <- model.matrix(~ tab.order2$stage + 0,PathExpMat2)
    colnames(design) <- c("AD","CH")
    fit <- lmFit(PathExpMat2, design)
    cont.matrix <- makeContrasts(AD-CH, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
    head(tT,10)
    hist(tT$adj.P.Val)
    tT.select <- tT[ tT$adj.P.Val < 0.05,]
    nrow(tT.select)
    plotDensities(PathExpMat2,legend = F)


    
#clustering highly differentially expressed pathways
    
    library(mclust)
    
    dif.mat <- PathExpMat2[rownames(tT.select),]
    BIC <- mclustBIC(t(dif.mat))
    plot(BIC)
    summary(BIC)
  
    
 

    heatmap(as.matrix(dif.mat))  
    
    d  <- dist( data.matrix(t(dif.mat) ), method = "euclidian")
    hc <- hclust(d,method = "complete")
    
    plot(hc,labels=colnames((dif.mat)),cex=0.1)
    

# Try to run some PCA down the road.
    
    seq.pca <- prcomp(t(PathExpMat2), center = TRUE, scale. = TRUE)
    plot(seq.pca,type = "l")
    
    # The plot tells us after the 4th PC we can't get much additional variance explained
    
    scoreTot <- seq.pca$x[,1:3]
    # Let's work on 5 PCs
    
    library(ggfortify)
    
    set.seed(1)
    
    # Does not seem to be affected by study ... Yoohoo?
    autoplot(seq.pca, data = tab.order2, colour = "outcome")
    
    library(pca3d)
    
# This should be good for now.
#Clustering children

tab.order2$stage == "child"
child.mat <- PathExpMat2[,tab.order2$stage == "child"]


library(mclust)

BIC <- mclustBIC((child.mat))
plot(BIC)
summary(BIC)


mod.path     <- Mclust(t(child.mat), G = 5, modelName = "VEI")
clust.order  <- unlist(tapply(1:ncol((child.mat)),
                             as.factor(mod.path$classification),I)
                      ,use.names = F)
child.mat    <- data.matrix(child.mat)[,clust.order]


library(RColorBrewer)
# Just some color setting for visualization

my_group    <- as.numeric(as.factor(mod.path$classification))
my_col1     <- brewer.pal(8, "Set1")[my_group]
my_col1     <- my_col1[clust.order]




# Distance based clustering 

d  <- dist( data.matrix(t(child.mat) ), method = "euclidian")
hc <- hclust(d,method = "complete")

plot(hc,labels=colnames(child.mat),cex=0.5)

# Three seems suitable so lets just do the rest of coloring and heatmap   
hclusters  <- cutree(hc, h=33)
row.order  <- unlist(tapply(1:nrow(t(child.mat)),
                            as.factor(hclusters),I),use.names = F)

child.mat   <- child.mat[row.order,]

my_group   <- as.numeric(as.factor(hclusters))
my_col1    <- brewer.pal(9, "Set1")[my_group]
my_row_col <- brewer.pal(9, "Set1")[my_group]
my_row_col <- my_row_col[row.order]

jpeg("images/Model-basedmRNACluster.jpg")
heatmap(as.matrix(child.mat))
dev.off()

clust1 <- names(hclusters[hclusters ==1])
clust2 <- names(hclusters[hclusters ==2])
clust3 <- names(hclusters[hclusters ==3])


# Do I need to normalize the pathways? anwer me quickly

# run1.list$arrayID
# run1.set   <- sapply(exp.mat.sets)
# set.lists    <- sapply(exp.mat.sets, nrow)
# grep(paste(datasets, collapse = "|"), names(set.lists),value = T)
# 
# paste(datasets, collapse = "|")
