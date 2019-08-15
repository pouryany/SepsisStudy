rm(list = ls())

library(pheatmap)
library(limma)

run1.exp.mat <- readRDS("DATA/SepsisGPL570.rds")
run1.exp.mat <- run1.exp.mat  %>%mutate_if(.,is.numeric,function(X){X*X}) 
#%>%mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

PathExpTab  <- readRDS("DATA/pathGeneTable.rds")
PathExpTab  <- left_join(PathExpTab,run1.exp.mat, by = c("SYMBOL" = "Gene"))
PathExpTab  <- PathExpTab %>% group_by(., Pathway) %>% 
               summarise_if(is.numeric,mean, na.rm=T)

pdf("images/SquaredRank/unCorrectedPathwayIntesities.pdf")
plotDensities(PathExpTab[,2:ncol(PathExpTab)],legend = F,
              main = "Raw Pathway Intensity per Sample")
dev.off()


samp.tab    <- readxl::read_xlsx("TrueSepsisDataUsedforFinalAnalysisMay2015.xlsx")
tab.order   <- colnames(PathExpTab)
tab.order   <- tab.order[-c(1)]
description <- samp.tab %>% filter(., arrayID %in% tab.order) %>% 
    select(.,arrayID,stage,study,outcome)

tab.order  <- left_join(as_tibble(tab.order),
                        description, by= c("value"= "arrayID"))

PathExpMat <- as.matrix(PathExpTab[,2:ncol(PathExpTab)])

rownames(PathExpMat) <- pull(PathExpTab[,1])
PathExpMat <- as.data.frame(PathExpMat)

plotMA(PathExpMat)
## The following is differential pathway analysis

            design <- model.matrix(~ tab.order$stage + 0,PathExpMat)
            colnames(design) <- c("AD","CH")
            fit <- lmFit(PathExpMat, design)
            cont.matrix <- makeContrasts(AD-CH, levels=design)
            fit2 <- contrasts.fit(fit, cont.matrix)
            fit2 <- eBayes(fit2, 0.01)
            tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
            
            
            head(tT,20)



###### Just redoing things with a Z-normalization

PathExpMat2 <- as.data.frame(PathExpMat)
PathExpMat2 <- apply(PathExpMat2, 2, function(X){ (X - mean(X))/sd(X)})
PathExpMat2 <- as.data.frame(PathExpMat2)


plotMA(PathExpMat2[,3:5])

            design <- model.matrix(~ tab.order$stage + 0,PathExpMat2)
            colnames(design) <- c("AD","CH")
            fit <- lmFit(PathExpMat2, design)
            cont.matrix <- makeContrasts(AD-CH, levels=design)
            fit2 <- contrasts.fit(fit, cont.matrix)
            fit2 <- eBayes(fit2, 0.01)
            tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
            head(tT,20)
            hist(tT$adj.P.Val)
            tT.select <- tT[ tT$adj.P.Val < 0.0001 & abs(tT$logFC) > 2,]
            nrow(tT.select)
            
            labs <- as.data.frame(tab.order[,2:4])
            rownames(labs) <- tab.order$value
            pdf("images/SquaredRank/PathwayDE570.pdf")
            pheatmap((PathExpMat2[rownames(tT.select),]),show_rownames = F,
                     annotation_col = labs)
            dev.off()
            
            #Subtypes exist in children. Not much in adults. 
            
            
            
            
pdf("images/SquaredRank/NormalizedPathwayIntesities.pdf")
plotDensities(PathExpMat2,legend = F,
              main = "Z-normalized Pathway Intensity per Sample")
dev.off()



library(ggfortify)

# PCA of normalized samples
GPL570.pca <- prcomp(t(PathExpMat2), center = TRUE, scale. = TRUE)
plot(GPL570.pca,type = "l")

set.seed(1)
# Does not seem to be affected by study ... Yoohoo?
pdf("images/SquaredRank/PCAnormalizedGPL570.pdf")
autoplot(GPL570.pca, data = tab.order, colour = "study",shape = "stage")
dev.off()


colnames(PathExpMat2) == tab.order$value

#Batch effect for children is not that much of a problem seemingly

pdf("images/SquaredRank/ChildsHeatmapAllpaths.pdf")
pheatmap((PathExpMat2[,tab.order$stage == "child"]),show_rownames = F)
dev.off()

child.mat  <- (PathExpMat2[,tab.order$stage == "child"])
high.paths <- (apply(child.mat,1, var))

pdf("images/SquaredRank/ChildPathwayVariance.pdf")
plotDensities(high.paths, main = "Pathway Variances in Children GPL570")
dev.off()

child.mat  <- child.mat[high.paths > 1.5, ]

# Heatmap below is amazing. Does no show batch effects in big subsets
labs <-  tab.order[tab.order$stage == "child",]
labs <- as.data.frame(labs[,3:4])
rownames(labs) <- tab.order[tab.order$stage == "child",]$value

pdf("images/SquaredRank/ChildrenSubtypesGPL570.pdf")
pheatmap(child.mat, annotation_col = labs,show_rownames = F )
dev.off()

aa <- pheatmap(child.mat, annotation_col = labs,show_rownames = F, cutree_rows = 4 )
write.csv(cutree(aa$tree_row, k = 4), "SquaredRankchildClusters.CSV")


### What's next is to find translation of the subgroups

child.pca <- prcomp(t(child.mat),center = TRUE, scale. = TRUE)

pdf("images/SquaredRank/PCAnormalizedChilds570.pdf")
autoplot(child.pca, data = tab.order[tab.order$stage == "child",], colour = "study")
dev.off()






### Looking at adults. There not much subclasses here. We'll see what to do.


        adult.mat  <- (PathExpMat2[,tab.order$stage == "adult"])
        high.paths <- (apply(adult.mat,1, var))
        
        pdf("images/SquaredRank/adultPathwayVariance.pdf")
        plotDensities(high.paths, main = "Pathway Variances in adultren GPL570")
        dev.off()
        
        adult.mat  <- adult.mat[high.paths > .1, ]
        
        # Heatmap below is amazing. Does no show batch effects in big subsets
        labs <-  tab.order[tab.order$stage == "adult",]
        labs <- as.data.frame(labs[,3:4])
        rownames(labs) <- tab.order[tab.order$stage == "adult",]$value
        
        pdf("images/SquaredRank/adultrenSubtypesGPL570.pdf")
        pheatmap(adult.mat, annotation_col = labs,show_rownames = F )
        dev.off()
        
        aa <- pheatmap(adult.mat, annotation_col = labs,show_rownames = F, cutree_rows = 4 )
        write.csv(cutree(aa$tree_row, k = 4), "SquaredRankadultClusters.CSV")


        
        
        

# Next is to see the pathway subtypes and see what they mean. 
# The subset of high-variance pathways in children can cluster effectively
# For 570 Not much of batch effect. Both for 571 there might be big differences
# In the 571 study there is a big time difference. Big. 
# Have to set up a regression type of classification next. 

# PCA of uncorrected samples
GPL570.pca <- prcomp(t(PathExpMat), center = TRUE, scale. = TRUE)
set.seed(1)
# Does not seem to be affected by study ... Yoohoo?
pdf("images/SquaredRank/PCAuncorrected570.pdf")
autoplot(GPL570.pca, data = tab.order, colour = "study",shape = "stage")
dev.off()









## Working on GPL571
    run2.exp.mat <- readRDS("DATA/SepsisGPL571.rds")
    run2.exp.mat <- run2.exp.mat  %>% mutate_if(.,is.numeric,rank) 
#%>%mutate_if(.,is.numeric,function(X){X*X})
# Edit ranking scales above!

    PathGeneTab      <- readRDS("DATA/pathGeneTable.rds")
    PathExpTab.571  <- left_join(PathGeneTab,run2.exp.mat, by = c("SYMBOL" = "Gene"))
    PathExpTab.571   <- PathExpTab.571 %>% group_by(., Pathway) %>% 
                        summarise_if(is.numeric,mean, na.rm=T)

    
    pdf("images/SquaredRank/unCorrectedPathwayIntesities571.pdf")
    plotDensities(PathExpTab.571[,2:ncol(PathExpTab.571)],legend = F)
    dev.off()
    
    PathExpMat.571 <- as.matrix(PathExpTab.571[,2:ncol(PathExpTab.571)])
    
    rownames(PathExpMat.571) <- pull(PathExpTab.571[,1])
    PathExpMat2.571 <- as.data.frame(PathExpMat.571)
    PathExpMat2.571 <- apply(PathExpMat2.571, 2, function(X){ (X - mean(X))/sd(X)})
    PathExpMat2.571 <- as.data.frame(PathExpMat2.571)
    
    
    pdf("images/SquaredRank/CorrectedPathwayIntesities571.pdf")
    plotDensities(PathExpMat2.571,legend = F)
    dev.off()
    
        #Check PCA for batch
            GPL571.pca      <- prcomp(t(PathExpMat2.571))
            tab.order.571   <- colnames(PathExpMat2.571)
            description     <- samp.tab %>% filter(., arrayID %in% tab.order.571) %>% 
                select(.,arrayID,stage,study,outcome)
            
            tab.order.571  <- left_join(as_tibble(tab.order.571), description, by= c("value"= "arrayID"))

            pdf("images/SquaredRank/PCAnormalizedGPL571.pdf")
            autoplot(GPL571.pca, data = tab.order.571, colour = "study")
            dev.off()
            # There is this batch effect here that I cannot do anything about?
        
          
  
## Putting both studies together          
            
            big.tab <- samp.tab %>% filter(.,platform %in% c("GPL570", "GPL571")) %>% 
                         select(.,arrayID,stage,outcome,study,gender)
            big.tab2<- as.data.frame(big.tab[,2:4])
            rownames(big.tab2) <- big.tab$arrayID
            rownames(big.tab2) == colnames(Big.mat)
            
            Big.mat <- cbind(PathExpMat2,PathExpMat2.571)
            Big.mat <- Big.mat[,rownames(big.tab2)]
            #big.polish <- medpolish(Big.mat)
            #Big.mat    <- big.polish$residuals
            big.pc  <-  prcomp(t(Big.mat[,rownames(big.tab2)]),center = T,scale. = T)
            autoplot(big.pc, data = big.tab2, colour = "study")
            
            vars <- apply(Big.mat,1,var)
            hist(vars)
            
            
            pheatmap(Big.mat[vars > 1,],show_rownames = F, annotation_col = big.tab2,
                     clustering_method = "ward.D2")
            
             
            # Looking at all adults
            big.tab2$stage == "adult" 
            big.adult <- Big.mat[,  big.tab2$stage == "adult"]

            
            
            vars <- apply(big.adult,1,var)
            plotDensities(vars)
            
            pheatmap(big.adult[vars > 0.01,],show_rownames = F,clustering_method = "ward.D2",
                     annotation_col = big.tab2[big.tab2$stage == "adult",] )
            
            
            big.polish <- medpolish(big.adult)
            big.adult    <- big.polish$residuals
            plotDensities(big.adult)
            
            
            
            big.pc  <-  prcomp(t(big.adult),center = T,scale. = T)
            autoplot(big.pc, data = big.tab2[big.tab2$stage == "adult",], colour = "study")
            ### This is interesting. Three distinct clusters. Adults are very similar?
            ### It's hard to get distinction by groups
            ### The batch effects disappear somewhat in pca but show up hierarchical
            ### Comparison of two plots says a lot
            ### Nothing Nothing. Batch dominates?!?!
            
# PCA Analysis too inspect the batch effects.

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
