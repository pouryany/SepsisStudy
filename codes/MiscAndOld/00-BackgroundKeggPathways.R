# The script to generate a table of pathway genes. This is based on the list in
# pathprint. Downloaded as of May 30th 2019. 


rm(list = ls())
load("DATA/pathprint.Hs.gs.rda")
load("DATA/genesets.rda")
 
pathList  <- pathprint.Hs.gs
pathList2 <- lapply(pathList, function(X){clusterProfiler::bitr(X,"ENTREZID",
                                                      "SYMBOL", 
                                                      OrgDb =  org.Hs.eg.db)} )


pathListTranslate <-list()

for(i in 1:length(pathList2)){
    temp <- data.frame(Pathway = names(pathList2[i]),pathList2[[i]])
    pathListTranslate <- rbind(pathListTranslate,temp)
    
    }
pathGeneTable <- as_tibble(pathListTranslate)
saveRDS(pathGeneTable,"DATA/pathGeneTable.rds")

