# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("topGO")
# BiocManager::install("Rgraphviz")


library(topGO)
library(ggplot2)
library(Rgraphviz)


setwd("E:/Splicata/16S-RNA-DNA/RNA_function")

# b) Lista de genes significativos
sigGenes <- read.table("list_corr_genes.txt")$V1
sigGenes <- unique(sigGenes)

# a) Lista de genes del universo
nocorr <- read.table("list_nocorr_genes.txt")$V1

allGenes <- c(sigGenes,nocorr)

# Vector binario nombre → 0/1
geneList <- factor(as.integer(allGenes %in% sigGenes))
names(geneList) <- allGenes

# c) Mapeo gene → GO
geneID2GO <- readMappings("gene2GO_comma_Splicata.map")

# d) Crear objeto topGO
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Tabla de resultados
allResFish <- GenTable(GOdata, 
                   classic = resultClassic, 
                   weight01 = resultWeight01,
                   elim = resultElim,
                   orderBy = "classic",
                   topNodes = 1256)

write.table(allResFish, "GO_corr_BP_Fisher.txt", sep="\t", quote=FALSE, row.names=FALSE)

allResFish <- subset(allResFish, Significant > Expected)
allResFish$statistic <- "Fisher"

claResFish <- allResFish[which(allResFish$classic < 0.05 & allResFish$Annotated > 3), ]
weiResFish <- allResFish[which(allResFish$weight01 < 0.05 & allResFish$Annotated > 3), ]
eliResFish <- allResFish[which(allResFish$elim < 0.05 & allResFish$Annotated > 3), ]
  
#   
#   
#   
# resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "KS")
# resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "KS")
# resultElim <- runTest(GOdata, algorithm = "elim", statistic = "KS")
# 
# allResKS <- GenTable(GOdata, 
#                        classic = resultClassic, 
#                        weight01 = resultWeight01,
#                        elim = resultElim,
#                        orderBy = "classic",
#                        topNodes = 1256)
# 
# write.table(allResKS, "GO_corr_BP_KS.txt", sep="\t", quote=FALSE, row.names=FALSE)
# 
# allResKS <- subset(allResKS, Significant > Expected)
# allResKS$statistic <- "KS" 
# 
# claResKS <- allResKS[which(allResKS$classic < 0.05 & allResKS$Annotated > 3), ]
# weiResKS <- allResKS[which(allResKS$weight01 < 0.05 & allResKS$Annotated > 3), ]
# eliResKS <- allResKS[which(allResKS$elim < 0.05 & allResKS$Annotated > 3), ]


claRes <- rbind(claResFish)
weiRes <- rbind(weiResFish)
eliRes <- rbind(eliResFish)



library(ggplot2)

pdf("classic_corr_BP.pdf",height=4, width=4.5)
ggplot(claRes, aes(x=reorder(Term, -as.numeric(classic)), y=as.numeric(classic))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(classic), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
    labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("weight01_corr_BP.pdf",height=4, width=4.5)
ggplot(weiRes, aes(x=reorder(Term, -as.numeric(weight01)), y=as.numeric(weight01))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(weight01), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("elim_corr_BP.pdf",height=4, width=4.5)
ggplot(eliRes, aes(x=reorder(Term, -as.numeric(elim)), y=as.numeric(elim))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(elim), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()



######

library(topGO)
library(ggplot2)
library(Rgraphviz)


setwd("E:/Splicata/16S-RNA-DNA/RNA_function")

# b) Lista de genes significativos
sigGenes <- read.table("list_corr_genes.txt")$V1
sigGenes <- unique(sigGenes)

# a) Lista de genes del universo
nocorr <- read.table("list_nocorr_genes.txt")$V1

allGenes <- c(sigGenes,nocorr)

# Vector binario nombre → 0/1
geneList <- factor(as.integer(allGenes %in% sigGenes))
names(geneList) <- allGenes

# c) Mapeo gene → GO
geneID2GO <- readMappings("gene2GO_comma_Splicata.map")

# d) Crear objeto topGO
GOdata <- new("topGOdata",
              ontology = "MF",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Tabla de resultados
allResFish <- GenTable(GOdata, 
                       classic = resultClassic, 
                       weight01 = resultWeight01,
                       elim = resultElim,
                       orderBy = "classic",
                       topNodes = 764)

write.table(allResFish, "GO_corr_MF_Fisher.txt", sep="\t", quote=FALSE, row.names=FALSE)

allResFish <- subset(allResFish, Significant > Expected)
allResFish$statistic <- "Fisher"

claResFish <- allResFish[which(allResFish$classic < 0.05 & allResFish$Annotated > 3), ]
weiResFish <- allResFish[which(allResFish$weight01 < 0.05 & allResFish$Annotated > 3), ]
eliResFish <- allResFish[which(allResFish$elim < 0.05 & allResFish$Annotated > 3), ]
# 
# resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "KS")
# resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "KS")
# resultElim <- runTest(GOdata, algorithm = "elim", statistic = "KS")
# 
# allResKS <- GenTable(GOdata, 
#                      classic = resultClassic, 
#                      weight01 = resultWeight01,
#                      elim = resultElim,
#                      orderBy = "classic",
#                      topNodes = 764)
# 
# write.table(allResKS, "GO_corr_MF_KS.txt", sep="\t", quote=FALSE, row.names=FALSE)
# 
# allResKS <- subset(allResKS, Significant > Expected)
# allResKS$statistic <- "KS" 
# 
# claResKS <- allResKS[which(allResKS$classic < 0.05 & allResKS$Annotated > 3), ]
# weiResKS <- allResKS[which(allResKS$weight01 < 0.05 & allResKS$Annotated > 3), ]
# eliResKS <- allResKS[which(allResKS$elim < 0.05 & allResKS$Annotated > 3), ]


claRes <- rbind(claResFish)
weiRes <- rbind(weiResFish)
eliRes <- rbind(eliResFish)


library(ggplot2)

pdf("classic_corr_MF.pdf",height=4, width=4.5)
ggplot(claRes, aes(x=reorder(Term, -as.numeric(classic)), y=as.numeric(classic))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(classic), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("weight01_corr_MF.pdf",height=4, width=4.5)
ggplot(weiRes, aes(x=reorder(Term, -as.numeric(weight01)), y=as.numeric(weight01))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(weight01), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("elim_corr_MF.pdf",height=4, width=4.5)
ggplot(eliRes, aes(x=reorder(Term, -as.numeric(elim)), y=as.numeric(elim))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(elim), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()





##### 


# b) Lista de genes significativos
sigGenes <- read.table("list_corr_genes_onlychr.txt")$V1

# a) Lista de genes del universo
nocorr <- read.table("list_nocorr_genes_onlychr.txt")$V1

allGenes <- c(sigGenes,nocorr)

# Vector binario nombre → 0/1
geneList <- factor(as.integer(allGenes %in% sigGenes))
names(geneList) <- allGenes

# c) Mapeo gene → GO
geneID2GO <- readMappings("gene2GO_comma_Splicata.map")

# d) Crear objeto topGO
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# Tabla de resultados
allResFish <- GenTable(GOdata, 
                       classic = resultClassic, 
                       weight01 = resultWeight01,
                       elim = resultElim,
                       orderBy = "classic",
                       topNodes = 1256)

write.table(allResFish, "GO_corr_onlychr_Fisher.txt", sep="\t", quote=FALSE, row.names=FALSE)

allResFish <- subset(allResFish, Significant > Expected)
allResFish$statistic <- "Fisher"

claResFish <- allResFish[which(allResFish$classic < 0.05 & allResFish$Annotated > 3), ]
weiResFish <- allResFish[which(allResFish$weight01 < 0.05 & allResFish$Annotated > 3), ]
eliResFish <- allResFish[which(allResFish$elim < 0.05 & allResFish$Annotated > 3), ]

# 
# 
# 
# resultWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "KS")
# resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "KS")
# resultElim <- runTest(GOdata, algorithm = "elim", statistic = "KS")
# 
# allResKS <- GenTable(GOdata, 
#                      classic = resultClassic, 
#                      weight01 = resultWeight01,
#                      elim = resultElim,
#                      orderBy = "classic",
#                      topNodes = 1245)
# 
# write.table(allResKS, "GO_corr_onlychr_KS.txt", sep="\t", quote=FALSE, row.names=FALSE)
# 
# allResKS <- subset(allResKS, Significant > Expected)
# allResKS$statistic <- "KS" 
# 
# claResKS <- allResKS[which(allResKS$classic < 0.05 & allResKS$Annotated > 3), ]
# weiResKS <- allResKS[which(allResKS$weight01 < 0.05 & allResKS$Annotated > 3), ]
# eliResKS <- allResKS[which(allResKS$elim < 0.05 & allResKS$Annotated > 3), ]


claRes <- rbind(claResFish)
weiRes <- rbind(weiResFish)
eliRes <- rbind(eliResFish)


library(ggplot2)

pdf("classic_corr_onlychr.pdf",height=4, width=4.5)
ggplot(claRes, aes(x=reorder(Term, -as.numeric(classic)), y=as.numeric(classic))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(classic), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("weight01_corr_onlychr.pdf",height=4, width=4.5)
ggplot(weiRes, aes(x=reorder(Term, -as.numeric(weight01)), y=as.numeric(weight01))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(weight01), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

pdf("elim_corr_onlychr.pdf",height=4, width=4.5)
ggplot(eliRes, aes(x=reorder(Term, -as.numeric(elim)), y=as.numeric(elim))) +
  geom_segment( lwd=1, aes(x=Term, xend=Term, y=0, yend=as.numeric(elim), color=Significant/Expected)) + 
  geom_point(aes(color = Significant/Expected, pch=statistic, size = Annotated)) + 
  coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(x="GO Term", y="p-value")+
  ylim(values=c(0,0.05))+
  scale_fill_gradient(low = "lightblue1", high = "navy") +
  scale_color_gradient(low = "lightblue1", high = "navy")
dev.off()

