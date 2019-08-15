rm(list =ls())
#######
##  GSE 28750 Sepsis data, whole blood: 10 Sepsis, 21 Normal,
##  11 Post surgical. Adults
######
 
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Wed May 29 14:54:35 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)

        #BiocManager::install("GEOquery")
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE28750", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- "00000000001111111111111111111122222222222"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G2-G0, G1-G0, G2-G1, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="none", number=Inf)
       

        ### This is the summarization part
        
        library(tidyverse)
        
        sum(is.na(tT$Gene.symbol))
        GSE287 <- exprs(gset)
        rownames(GSE287) <- tT$Gene.symbol
        GSE287 <- as_tibble(GSE287, rownames = "Gene")
        
  
        GSE287 <- GSE287 %>% group_by(., Gene) %>% summarise_all(.,mean)
        GSE287 <- GSE287 %>% filter (., Gene != "")
        
        #Table having samples in column, pathway-genes in rows
        PathGeneTab <- readRDS("DATA/pathGeneTable.rds")
        pathGeneExp <- left_join(PathGeneTab,GSE287, by = c("SYMBOL" = "Gene"))
        
        PathExpTab  <- pathGeneExp %>% group_by(., Pathway) %>% 
                       summarise_if(is.numeric,mean, na.rm=T)
        

    # Some plotting. Pathways seems pretty normalized.
        
        plot(density(pull(PathExpTab[,5])))
        
        plotDensities(gset, legend = F)
        plotMA(gset)
        
    # Some pathway level expression DE analysis. 
    # Should move to non-parametric.
        PathExpMat           <- as.matrix.data.frame(PathExpTab[,2:42])
        rownames(PathExpMat) <- as.character(pull(PathExpTab[,1]))
        fit <- lmFit(PathExpMat, design)
        cont.matrix <- makeContrasts(G2-G0, G1-G0, G2-G1, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="F", number=Inf)
        
        head(tT)
        plot(density(tT$P.Value))
        
    # Ok, there seems to be missing ones. do I remove them simply? or do I 
    # match the lists across all of the assays/pathways.       
        grep("Pentos",pathGeneExp$Pathway,value = T)
        test.this <- pathGeneExp %>%
                     filter(.,Pathway == "Pentose phosphate pathway (KEGG)")

        test.this[is.na(test.this$GSM712479),]
        
        
 
        
#######
##  GSE 13015 Sepsis data, whole blood: 29 Sepsis, 10 Normal
##  Adults. The data have treatments/agents. It also has pathogens
######

        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Wed May 29 15:03:49 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE13015", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- "000000000000000111111111100000000000000"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
        write.table(tT, file=stdout(), row.names=F, sep="\t")
        
        
        ################################################################
        #   Boxplot for selected GEO samples
        library(Biobase)
        library(GEOquery)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE13015", GSEMatrix =TRUE, getGPL=FALSE)
        if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # group names for all samples in a series
        gsms <- "000000000000000111111111100000000000000"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        sml <- paste("G", sml, sep="") # set group names
        
        # order samples by group
        ex <- exprs(gset)[ , order(sml)]
        sml <- sml[order(sml)]
        fl <- as.factor(sml)
        labels <- c("Sepsis","Normal","PostSurgical")
        

        
#######
##  GSE 10474 Sepsis data, whole blood: 13 Lung Injury Sepsis Sepsis, 21 Sepsis
##  Adults. 
######        
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Wed May 29 15:15:11 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE10474", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- "3333333333333000000000000000000000"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G3-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
        
#######
##  GSE 40586 Sepsis data, whole blood: 14 Sepsis, 18 Normal
##  Adults. Check the difference in the number of sample with the paper.
        ## Group three is sepsis?! Are there children in here?
######    

        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Wed May 29 15:29:06 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE40586", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- "3003X00X30003030XX30X3X3X30300303303000"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # eliminate samples marked as "X"
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[ ,sel]
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G3-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
       

#######
##  GSE57065 Sepsis data, whole blood: 28 Sepsis at hour 0, 25 non sepsis
##  There are other sepsis samples, 
##  Adults. There is male/female variable as well.
######           

        
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:11:13 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE57065", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- paste0("1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1XX1X",
                       "X1X1XX1XX1XX1XX1X1XX1XX1XX1XX1XX000000000000000000",
                       "0000000")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # eliminate samples marked as "X"
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[ ,sel]
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        
#######
## GSE33341 Sepsis data, Adults, whole blood: 51 Sepsis, 43 Normal
##  There are processing batch variables. as well as agent/race/gender/infectio site
######        
        
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:18:09 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE33341", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- paste0("11111111111111111111111111111111111111110000000000",
                       "00000000000000000000000000000000011111111111")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

#
#
#
#######
## GSE4607 Sepsis data, Confusing subsets. Should ask soon.
######                  
#
#
#

        
                
        
#
#
#
#######
## GSE9692 Sepsis data, Children. The number of samples in the data and
## the publication are different. Ask why?
######                  
#
#
#
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:34:40 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE9692", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- "111111111111111111111000000000001111111100001"
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        

#
#
#
#######
## GSE26440 Sepsis data, Children. Sepsis 98, Normal 32 The number of samples in the data and
## the publication are different. Ask why?
######                  
#
#
#
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:39:13 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE26440", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- paste0("11111111111111111111111111111111111111110001111000",
                       "00111111111111111111111111110010000100001011110000",
                       "000001111000101111111111111111")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        

#
#
#
#######
## GSE26378 Sepsis data, Children. Sepsis 82, Normal 21 The number of samples in the data and
## the publication are different. Ask why?
######                  
#
#
#
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:53:41 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE26378", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- paste0("11111100100001111111111111111111111111101111111101",
                       "11111111111100111111111111100000111111111111000000",
                       "111")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        plotDensities(gset, legend = F)
        pData(gset)
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G1-G0, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        
        


#
#
#
#######
## GSE139 Sepsis data, Children. Samples: Normal 18, Sepsis 52, Septic Shock 106, SIRS 27
######                  
#
#
#

        
        
        
        # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
        # R scripts generated  Thu May 30 13:58:22 EDT 2019
        
        ################################################################
        #   Differential expression analysis with limma
        library(Biobase)
        library(GEOquery)
        library(limma)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE13904", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group names for all samples
        gsms <- paste0("000000000000000000333333333333333333333333333XXXXX",
                       "XXXXXXXXXXXXXXXXXXX1111111111111111111111111111111",
                       "11111111111111111111122222222222222222222222222222",
                       "22222222222222222222222222222222222222222222222222",
                       "222222222222222222222222222")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        
        # eliminate samples marked as "X"
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[ ,sel]
        
        # log2 transform
        ex <- exprs(gset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset) <- log2(ex) }
        
        # set up the data and proceed with analysis
        sml <- paste("G", sml, sep="")    # set group names
        fl <- as.factor(sml)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(G3-G0, G1-G0, G2-G1, G3-G2, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
        write.table(tT, file=stdout(), row.names=F, sep="\t")
        
        
        ################################################################
        #   Boxplot for selected GEO samples
        library(Biobase)
        library(GEOquery)
        
        # load series and platform data from GEO
        
        gset <- getGEO("GSE13904", GSEMatrix =TRUE, getGPL=FALSE)
        if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # group names for all samples in a series
        gsms <- paste0("000000000000000000333333333333333333333333333XXXXX",
                       "XXXXXXXXXXXXXXXXXXX1111111111111111111111111111111",
                       "11111111111111111111122222222222222222222222222222",
                       "22222222222222222222222222222222222222222222222222",
                       "222222222222222222222222222")
        sml <- c()
        for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
        sml <- paste("G", sml, sep="")  #set group names
        
        # eliminate samples marked as "X"
        sel <- which(sml != "X")
        sml <- sml[sel]
        gset <- gset[ ,sel]
        
        # order samples by group
        ex <- exprs(gset)[ , order(sml)]
        sml <- sml[order(sml)]
        fl <- as.factor(sml)
        labels <- c("Normal","Sepsis","Septic+Shock","SIRS")
        
        