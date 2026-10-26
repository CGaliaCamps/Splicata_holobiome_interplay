####load your packages and set directory####

# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("DESeq2", force=TRUE)
# BiocManager::install("apeglm", force=TRUE)
# BiocManager::install("vsn", force=TRUE)
# install.packages("stringr")
# install.packages("svglite")


library(ggplot2)
library(factoextra)

library(apeglm)
library(vsn)
library(stringr)

library(svglite)
library(pheatmap)
library(clusterSim)

library(DESeq2)

setwd("E:/Splicata/16S-RNA-DNA/1Dif_express/New_stuff_15")

####      load your matrix and metadata####
cts <- as.data.frame(read.table("Splicata_counts_15.txt", header=T)) 
colnames(cts) <- sub("_15_sorted.bam", "", colnames(cts))
cts <- cts[rowSums(cts)>0,] 
cts <- as.matrix(cts)
names <- rownames(cts)
#necessitem que els valors de la matriu estiguin amb format integer
class(cts[1,]) 
cts <- apply(cts, MARGIN = 2, FUN = as.integer) 
class(cts[1,]) 
rownames(cts) <- names#aquesta funció és perque al transformar perdem el nom de les files, aixi els recuperem copiant-los de la matriu anterior
#comprovar abans si l'ordre dels gens és el mateix i no sha modificat

#carreguem la matriu amb les metadades, la informació dels individus
coldata <- read.table(file="metadata.txt", header=T)
rownames(coldata) <- coldata$Sample
coldata <- coldata[,2:6]
coldata[,1:5] <-lapply(coldata[,1:5], factor)
str(coldata)

#ordenem els individus de forma alfabetica per tenir les dues matrius en el mateix ordre
coldata <- coldata[order(rownames(coldata)),c("Pop","Chr2", "Chr4","Chr11", "Chr16")]

#ordenem
cts <- cts[, rownames(coldata)] #Ordena les dues matrius igual
all(rownames(coldata) == colnames(cts)) #aqui indica que estan e el mateix ordre



#generate annotation files for microbiome (necessary in next steps)
gff<-read.table("splicata_round3_lnc_sno_trna_notracks.gtf", sep="\t")
gff <- gff[which(gff$V3=="transcript"),]
gff <- gff[,c(1,4,5,9)]

colnames(gff) <- c("Chrom","start","end", "gene")

GOTERM <- gff #Next in the script

gff$gene <- gsub("; .*", "", gff$gene)
gff$gene <- gsub("transcript_id ", "", gff$gene)
gff$gene <- gsub("\\.[^.]*$", "", gff$gene)
gff$gene <- gsub("-.*", "", gff$gene)

gff <- gff[grepl("^Chrom*", gff$Chrom) ,]

gff <- gff[!duplicated(gff),]

####      Do a PCA of your samples to check possible weird samples####
cts.normal <- as.data.frame(t(cts))
cts.normal <- data.Normalization(cts.normal, type = "n12", normalization = "column")
RNA_PCA <- prcomp(cts.normal)
summary(RNA_PCA)

write.table(cts.normal, file="Splicata_counts_15_normalized.txt", sep="\t", row.names=TRUE) 

RNA_PCA_plot <- as.data.frame(RNA_PCA$x)
RNA_PCA_plot <- cbind(RNA_PCA_plot, coldata)
PCA_loadings <- as.data.frame(unclass(RNA_PCA$rotation))
eig.val <- get_eigenvalue(RNA_PCA)
#screeplot(RNA_PCA, type = "barplot")

pdf("PCA_alldata_RNA.pdf", height=4.55, width=5)
ggplot(data = RNA_PCA_plot, aes(x = PC1, y = PC2, color=Pop))+
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=0.75)+
  geom_point(data = RNA_PCA_plot, size = 5, pch= 19,stroke = 1.5, alpha=0.8)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
    scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(eig.val[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(eig.val[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

####


####      Differential expression models######

ddsPop <- DESeqDataSetFromMatrix(countData = cts,
                                 colData = coldata,
                                 design = ~ Pop) #without interactions as we do not have all levels interacting in most of the cases
ddsChr2 <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ Chr2)
ddsChr4 <- DESeqDataSetFromMatrix(countData = cts[,-14],
                                  colData = coldata[-14,],
                                  design = ~ Chr4)
ddsChr11 <- DESeqDataSetFromMatrix(countData = cts[,-c(8,10)],
                                   colData = coldata[-c(8,10),],
                                   design = ~ Chr11)
ddsChr16 <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design = ~ Chr16)

#####
#####
#####
####      Differential expression (populations and inversions) using all genes######


#set the references for each factor
ddsPop$Pop <- relevel(ddsPop$Pop, ref = "BLA")
ddsChr2$Chr2 <- relevel(ddsChr2$Chr2, ref = "AB")
ddsChr4$Chr4 <- relevel(ddsChr4$Chr4, ref = "AB")
ddsChr11$Chr11 <- relevel(ddsChr11$Chr11, ref = "AB")
ddsChr16$Chr16 <- relevel(ddsChr16$Chr16, ref = "AB")

 genesToRemove <- which(!rownames(ddsPop) %in% gff$gene)
 
 ddsPop<-ddsPop[-genesToRemove,]
 ddsChr2<-ddsChr2[-genesToRemove,]
 ddsChr4<-ddsChr4[-genesToRemove,]
 ddsChr11<-ddsChr11[-genesToRemove,]
 ddsChr16<-ddsChr16[-genesToRemove,]

#analyze differential expression!
ddseqPop <- DESeq(ddsPop)
ddseqChr2 <- DESeq(ddsChr2)
ddseqChr4 <- DESeq(ddsChr4)
ddseqChr11 <- DESeq(ddsChr11)
ddseqChr16 <- DESeq(ddsChr16)

#set p-value as 0.05 to be significantly different and order according to it
resPop_95 <- results(ddseqPop, alpha=0.05)
resChr2_95 <- results(ddseqChr2, alpha=0.05)
resChr4_95 <- results(ddseqChr4, alpha=0.05)
resChr11_95 <- results(ddseqChr11, alpha=0.05)
resChr16_95 <- results(ddseqChr16, alpha=0.05)

resPop_95 <- resPop_95[order(resPop_95$pvalue),]
resChr2_95 <- resChr2_95[order(resChr2_95$pvalue),]
resChr4_95 <- resChr4_95[order(resChr4_95$pvalue),]
resChr11_95 <- resChr11_95[order(resChr11_95$pvalue),]
resChr16_95 <- resChr16_95[order(resChr16_95$pvalue),]

#get summaries and export tables
summary(resPop_95)
summary(resChr2_95)
summary(resChr4_95)
summary(resChr11_95)
summary(resChr16_95)

write.table(resPop_95, file="DESeq_Pop_95.txt", sep="\t", row.names=TRUE)
write.table(resChr2_95, file="DESeq_Chr2_95.txt", sep="\t", row.names=TRUE)
write.table(resChr4_95, file="DESeq_Chr4_95.txt", sep="\t", row.names=TRUE)
write.table(resChr11_95, file="DESeq_Chr11_95.txt", sep="\t", row.names=TRUE)
write.table(resChr16_95, file="DESeq_chr16_95.txt", sep="\t", row.names=TRUE)
####      Identify genes expressed differently and represent them ####

#p adjust value below 0.05 and with a mean number of reads of 2
Sig_Pop_95 <- subset(resPop_95, padj < 0.05 & baseMean > 2) ## NO DIFFERENCES BETWEEN baseMean 2 & 10
Sig_Chr2_95 <- subset(resChr2_95, padj < 0.05 & baseMean > 2 )
Sig_Chr4_95 <- subset(resChr4_95, padj < 0.05 & baseMean > 2 )
Sig_Chr11_95 <- subset(resChr11_95, padj < 0.05 & baseMean > 2 )
Sig_Chr16_95 <- subset(resChr16_95, padj < 0.05 & baseMean > 2 )

summary(Sig_Pop_95)
summary(Sig_Chr2_95)
summary(Sig_Chr4_95)
summary(Sig_Chr11_95)
summary(Sig_Chr16_95)

Sel_Pop <- rownames(Sig_Pop_95)
Sel_Chr2 <-rownames(Sig_Chr2_95)
Sel_Chr4 <-rownames(Sig_Chr4_95)
Sel_Chr11 <-rownames(Sig_Chr11_95)
Sel_Chr16 <-rownames(Sig_Chr16_95)

#get an upset plot
list_all <- list(Pop= Sel_Pop, 
                 Chr2= Sel_Chr2, 
                 Chr4= Sel_Chr4,
                 Chr11=Sel_Chr11,
                 Chr16=Sel_Chr16
                 )

library(ggVennDiagram)
venn_RNAexp<- ggVennDiagram(list_all, label= "c", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none")

ggsave(venn_RNAexp, file="Venn_RNAexp.svg", device = "svg", height=5, width=5)

names<-names(list_all)

list_all <- t(splitstackshape:::charMat(list_all, fill = 0, mode="binary"))
colnames(list_all) <- names 
list_all <- rev(as.data.frame(list_all))
sets <-    colnames(list_all) 

pdf("RNA_upset.pdf", width=6, height=5.5)
UpSetR::upset(list_all, sets = sets, mb.ratio = c(0.7, 0.3), 
             # sets.bar.color = c("purple","#b05446","#9dcfdd"), 
              order.by = c("freq","degree"), keep.order=T, line.size = 0.8
              ,sets.x.label = "Number of Diff. Expressed RNA", mainbar.y.label = "Number of Diff. Expressed RNA")
dev.off()



#Remove the genes in common (interactions?)
vector1 <- intersect(Sel_Pop, Sel_Chr2)
vector2 <- intersect(Sel_Pop, Sel_Chr4)
vector3 <- intersect(Sel_Pop, Sel_Chr11)
vector4 <- intersect(Sel_Pop, Sel_Chr16)
vector5 <- intersect(Sel_Chr2, Sel_Chr4)
vector6 <- intersect(Sel_Chr2, Sel_Chr11)
vector7 <- intersect(Sel_Chr2, Sel_Chr16)
vector8 <- intersect(Sel_Chr4, Sel_Chr11)
vector9 <- intersect(Sel_Chr4, Sel_Chr16)
vector10 <- intersect(Sel_Chr11, Sel_Chr16)

interact <- c(vector1, vector2, vector3, vector4,vector5,vector6,vector7,vector8,vector9,vector10)
interact <- unique(interact)

#We keep only those genes we know are specific for the factor we tested
Sel_Pop <- setdiff(Sel_Pop, interact)
Sel_Chr2 <- setdiff(Sel_Chr2, interact)
Sel_Chr4 <- setdiff(Sel_Chr4, interact)
Sel_Chr11 <- setdiff(Sel_Chr11, interact)
Sel_Chr16 <- setdiff(Sel_Chr16, interact)

#Export the gene names that are differentially expressed
write.table(Sel_Pop, file="diff_pop_genes.txt", row.names=F, col.names=F)
write.table(Sel_Chr2, file="diff_chr2_genes.txt", row.names=F, col.names=F)
write.table(Sel_Chr4, file="diff_chr4_genes.txt", row.names=F, col.names=F)
write.table(Sel_Chr11, file="diff_chr11_genes.txt", row.names=F, col.names=F)
write.table(Sel_Chr16, file="diff_chr16_genes.txt", row.names=F, col.names=F)

####      Plot heatmaps####

#normalize our read matrix 
norm_Pop <- normTransform(ddseqPop)
norm_Chr2 <- normTransform(ddseqChr2)
norm_Chr4 <- normTransform(ddseqChr4)
norm_Chr11 <- normTransform(ddseqChr11)
norm_Chr16 <- normTransform(ddseqChr16)

nomral_counts<-assay(norm_Pop)
write.table(nomral_counts, file="normal_counts.txt", sep="\t", row.names=TRUE)

#generate annotation files for the pheatmap and give them proper rownames
coldata4 <- coldata[-14,]
coldata11 <- coldata[-c(8,10),]

annotation_Pop <- data.frame(Pop=coldata$Pop)
annotation_Chr2 <- data.frame(Chr2=coldata$Chr2, Pop=coldata$Pop)
annotation_Chr4 <- data.frame(Chr4=coldata4$Chr4, Pop=coldata4$Pop)
annotation_Chr11 <- data.frame(Chr11=coldata11$Chr11, Pop=coldata11$Pop)
annotation_Chr16 <- data.frame(Chr16=coldata$Chr16, Pop=coldata$Pop)

rownames(annotation_Pop) <- rownames(annotation_Chr2) <-  rownames(annotation_Chr16) <- rownames(coldata)
rownames(annotation_Chr4) <- rownames(coldata4) 
rownames(annotation_Chr11) <- rownames(coldata11) 

#generate color file for our annotation
colors_Pop = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                  Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue",Chromosome_5="black",
                            Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                            Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                            Chromosome_16="lightblue"))
colors_Chr2 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr2 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                             Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                             Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                             Chromosome_16="lightblue4"))
colors_Chr4 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr4 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue",Chromosome_5="black",
                             Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                             Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                             Chromosome_16="lightblue4"))
colors_Chr11 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr11 = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
                    Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                              Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                              Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                              Chromosome_16="lightblue4"))
colors_Chr16 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr16 = c(AA="#ecd6c0", AB="#834f23"),
                    Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                              Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                              Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                              Chromosome_16="lightblue"))


gff_map <- unique(gff[, c("gene", "Chrom")])
rownames(gff_map) <- gff_map$gene
annotation_row_pop <- data.frame(
  Chrom = gff_map[Sel_Pop, "Chrom"],
  row.names = Sel_Pop,
  stringsAsFactors = FALSE
)

annotation_row_chr2 <- data.frame(
  Chrom = gff_map[Sel_Chr2, "Chrom"],
  row.names = Sel_Chr2,
  stringsAsFactors = FALSE
)

annotation_row_chr4 <- data.frame(
  Chrom = gff_map[Sel_Chr4, "Chrom"],
  row.names = Sel_Chr4,
  stringsAsFactors = FALSE
)

annotation_row_chr11 <- data.frame(
  Chrom = gff_map[Sel_Chr11, "Chrom"],
  row.names = Sel_Chr11,
  stringsAsFactors = FALSE
)

annotation_row_chr16 <- data.frame(
  Chrom = gff_map[Sel_Chr16, "Chrom"],
  row.names = Sel_Chr16,
  stringsAsFactors = FALSE
)






#PLOT. JUST DO IT!!!!


callback_pop = function(hc, mat){
  sv = svd(t(mat))$v[,2]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Pop.pdf", width=7, height =6)
pheatmap(assay(norm_Pop)[Sel_Pop,],
         cluster_rows=T, show_rownames=F, cluster_cols=T, scale="row",
         annotation_col=annotation_Pop, annotation_colors = colors_Pop, annotation_row = annotation_row_pop,
         clustering_callback=callback_pop)
dev.off()


callback_chr2 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr2.pdf", width=7, height =5)
pheatmap(assay(norm_Chr2)[Sel_Chr2,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr2, annotation_colors = colors_Chr2, annotation_row = annotation_row_chr2,
         clustering_callback=callback_chr2)
dev.off()


callback_chr4 = function(hc, mat){
  sv = svd(t(mat))$v[,14]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr4.pdf", width=7, height =5)
pheatmap(assay(norm_Chr4)[Sel_Chr4,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr4, annotation_colors = colors_Chr4, annotation_row = annotation_row_chr4,
         )
dev.off()


callback_chr11 = function(hc, mat){
  sv = svd(t(mat))$v[,4]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr11.pdf", width=6.5, height =5)
pheatmap(assay(norm_Chr11)[Sel_Chr11,], 
         cluster_rows=T, show_rownames=F, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr11, annotation_colors = colors_Chr11, annotation_row = annotation_row_chr11,
         clustering_callback=callback_chr11)
dev.off()


callback_chr16 = function(hc, mat){
  sv = svd(t(mat))$v[,2]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr16_.pdf", width=7, height =5)
pheatmap(assay(norm_Chr16)[Sel_Chr16,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr16, annotation_colors = colors_Chr16, annotation_row = annotation_row_chr16,
         clustering_callback=callback_chr16)
dev.off()
dev.off()

####      Generate PCAs only with data from relevant differentially expessed genes####

#extract the matrices
PCA_pop <- t(assay(norm_Pop)[Sel_Pop,])
PCA_Chr2 <- t(assay(norm_Chr2)[Sel_Chr2,])
PCA_Chr4 <- t(assay(norm_Chr4)[Sel_Chr4,])
PCA_Chr11 <- t(assay(norm_Chr11)[Sel_Chr11,])
PCA_Chr16 <- t(assay(norm_Chr16)[Sel_Chr16,])

#do the principal component analyses
PCA_pop <- prcomp(PCA_pop)
PCA_Chr2 <- prcomp(PCA_Chr2)
PCA_Chr4 <- prcomp(PCA_Chr4)
PCA_Chr11 <- prcomp(PCA_Chr11)
PCA_Chr16 <- prcomp(PCA_Chr16)

#get your gene loadings to the PCAs
Pop_loadings <- as.data.frame(unclass(PCA_pop$rotation))
Chr2_loadings <- as.data.frame(unclass(PCA_Chr2$rotation))
Chr4_loadings <- as.data.frame(unclass(PCA_Chr4$rotation))
Chr11_loadings <- as.data.frame(unclass(PCA_Chr11$rotation))
Chr16_loadings <- as.data.frame(unclass(PCA_Chr16$rotation))

#Pop_loadings %>% arrange(desc(abs(PC1)))
#Get the % expained by each axis
Pop.eig <- get_eigenvalue(PCA_pop)
Chr2.eig <- get_eigenvalue(PCA_Chr2)
Chr4.eig <- get_eigenvalue(PCA_Chr4)
Chr11.eig <- get_eigenvalue(PCA_Chr11)
Chr16.eig <- get_eigenvalue(PCA_Chr16)

#generate a matrix to plot
PCA_pop <- as.data.frame(PCA_pop$x)
PCA_Chr2 <- as.data.frame(PCA_Chr2$x)
PCA_Chr4 <- as.data.frame(PCA_Chr4$x)
PCA_Chr11 <- as.data.frame(PCA_Chr11$x)
PCA_Chr16 <- as.data.frame(PCA_Chr16$x)

PCA_pop <- cbind(PCA_pop, coldata)
PCA_Chr2 <- cbind(PCA_Chr2, coldata)
PCA_Chr4 <- cbind(PCA_Chr4, coldata4)
PCA_Chr11 <- cbind(PCA_Chr11, coldata11)
PCA_Chr16 <- cbind(PCA_Chr16, coldata)


pdf("PCA_Pop_genes.pdf")
ggplot(data = PCA_pop, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Pop_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_pop, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
     # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
   #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
   labs(x = paste("PC1 (", round(Pop.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Pop.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr2_genes.pdf")
ggplot(data = PCA_Chr2, aes(x = PC1, y = PC2, color=Chr2))+
  geom_segment(data = Chr2_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr2), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr2, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23","#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23","#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr2.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr2.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr4_genes.pdf")
ggplot(data = PCA_Chr4, aes(x = PC1, y = PC2, color=Chr4))+
  geom_segment(data = Chr4_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr4), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr4, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr4.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr4.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr11_genes.pdf")
ggplot(data = PCA_Chr11, aes(x = PC1, y = PC2, color=Chr11))+
  geom_segment(data = Chr11_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(color=Chr11, fill=Chr11), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr11, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#834f23",  "#45322E"))+
  scale_fill_manual(values=c( "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr11.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr11.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr16_genes.pdf")
ggplot(data = PCA_Chr16, aes(x = PC1, y = PC2, color=Chr16))+
  geom_segment(data = Chr16_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr16), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr16, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr16.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr16.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()




#####
#####
#####
####      Differential expression (inversions) using chrom inversion genes######
ddseqChr2 <- DESeq(ddsChr2)
ddseqChr4 <- DESeq(ddsChr4)
ddseqChr11 <- DESeq(ddsChr11)
ddseqChr16 <- DESeq(ddsChr16)

gene_chr2 <- gff[which(gff$Chrom=="Chromosome_2" & gff$start>=1 & gff$end<=19480816),]
gene_chr4 <- gff[which(gff$Chrom=="Chromosome_4" & gff$start>=1 & gff$end<=16999389),]
gene_chr11 <- gff[which(gff$Chrom=="Chromosome_11" & gff$start>=1 & gff$end<=10440662),]
gene_chr16 <- gff[which(gff$Chrom=="Chromosome_16" & gff$start>=4371327 & gff$end<=18293200),]

#correct qvalues for only those genes in a  certain chromosome
#install.packages("FDRestimation")
library(FDRestimation)

pval_chr2 <- as.data.frame(results(ddseqChr2))
pval_chr2 <- pval_chr2[rownames(pval_chr2) %in% gene_chr2$gene, ]
pval_chr2 <- pval_chr2[!is.na(pval_chr2$pvalue), ]
pval_chr2_test <- pval_chr2$pvalue
FDR_Chr2 <-p.fdr(pvalues=pval_chr2_test, just.fdr=T)
pval_chr2$FDR <- FDR_Chr2
Sel_Chr2 <- rownames(subset(pval_chr2, FDR < 0.05 & baseMean > 2))

pval_chr4 <- as.data.frame(results(ddseqChr4))
pval_chr4 <- pval_chr4[rownames(pval_chr4) %in% gene_chr4$gene, ]
pval_chr4 <- pval_chr4[!is.na(pval_chr4$pvalue), ]
pval_chr4_test <- pval_chr4$pvalue
FDR_chr4 <-p.fdr(pvalues=pval_chr4_test, just.fdr=T)
pval_chr4$FDR <- FDR_chr4
Sel_Chr4 <- rownames(subset(pval_chr4, FDR < 0.05 & baseMean > 2))

pval_chr11 <- as.data.frame(results(ddseqChr11))
pval_chr11 <- pval_chr11[rownames(pval_chr11) %in% gene_chr11$gene, ]
pval_chr11 <- pval_chr11[!is.na(pval_chr11$pvalue), ]
pval_chr11_test <- pval_chr11$pvalue
FDR_chr11 <-p.fdr(pvalues=pval_chr11_test, just.fdr=T)
pval_chr11$FDR <- FDR_chr11
Sel_Chr11 <- rownames(subset(pval_chr11, FDR < 0.05 & baseMean > 2))

pval_chr16 <- as.data.frame(results(ddseqChr16))
pval_chr16 <- pval_chr16[rownames(pval_chr16) %in% gene_chr16$gene, ]
pval_chr16 <- pval_chr16[!is.na(pval_chr16$pvalue), ]
pval_chr16_test <- pval_chr16$pvalue
FDR_chr16 <-p.fdr(pvalues=pval_chr16_test, just.fdr=T)
pval_chr16$FDR <- FDR_chr16
Sel_Chr16 <- rownames(subset(pval_chr16, FDR < 0.05 & baseMean > 2))

#get an upset plot
list_all <- list(Chr2= Sel_Chr2, 
                 Chr4= Sel_Chr4,
                 Chr11=Sel_Chr11,
                 Chr16=Sel_Chr16
)

library(ggVennDiagram)
venn_RNAexp<- ggVennDiagram(list_all, label= "c", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none")

ggsave(venn_RNAexp, file="Venn_RNAexp_onlyChr.svg", device = "svg", height=5, width=5)

names<-names(list_all)

list_all <- t(splitstackshape:::charMat(list_all, fill = 0, mode="binary"))
colnames(list_all) <- names 
list_all <- rev(as.data.frame(list_all))
sets <-    colnames(list_all) 

pdf("RNA_upset_onlyChr.pdf", width=6, height=5.5)
UpSetR::upset(list_all, sets = sets, mb.ratio = c(0.7, 0.3), 
              # sets.bar.color = c("purple","#b05446","#9dcfdd"), 
              order.by = c("freq","degree"), keep.order=T, line.size = 0.8
              ,sets.x.label = "Number of Diff. Expressed RNA", mainbar.y.label = "Number of Diff. Expressed RNA")
dev.off()
dev.off()

#Remove the genes in common (interactions?)
vector5 <- intersect(Sel_Chr2, Sel_Chr4)
vector6 <- intersect(Sel_Chr2, Sel_Chr11)
vector7 <- intersect(Sel_Chr2, Sel_Chr16)
vector8 <- intersect(Sel_Chr4, Sel_Chr11)
vector9 <- intersect(Sel_Chr4, Sel_Chr16)
vector10 <- intersect(Sel_Chr11, Sel_Chr16)

interact <- c(vector5,vector6,vector7,vector8,vector9,vector10)
interact <- unique(interact)

#We keep only those genes we know are specific for the factor we tested
Sel_Chr2 <- setdiff(Sel_Chr2, interact)
Sel_Chr4 <- setdiff(Sel_Chr4, interact)
Sel_Chr11 <- setdiff(Sel_Chr11, interact)
Sel_Chr16 <- setdiff(Sel_Chr16, interact)

#Export the gene names that are differentially expressed
write.table(Sel_Chr2, file="diff_chr2_genes_onlychr.txt", row.names=F, col.names=F)
write.table(Sel_Chr4, file="diff_chr4_genes_onlychr.txt", row.names=F, col.names=F)
write.table(Sel_Chr11, file="diff_chr11_genes_onlychr.txt", row.names=F, col.names=F)
write.table(Sel_Chr16, file="diff_chr16_genes_onlychr.txt", row.names=F, col.names=F)

####      Plot heatmaps####

#normalize our read matrix 
norm_Chr2 <- normTransform(ddseqChr2)
norm_Chr4 <- normTransform(ddseqChr4)
norm_Chr11 <- normTransform(ddseqChr11)
norm_Chr16 <- normTransform(ddseqChr16)


colors_Chr2 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr2 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_2="lightblue"))
colors_Chr4 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr4 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_4="lightblue"))
colors_Chr11 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr11 = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
                    Chrom = c(Chromosome_11="lightblue"))
colors_Chr16 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr16 = c(AA="#ecd6c0", AB="#834f23"),
                    Chrom = c(Chromosome_16="lightblue"))


callback_chr2 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr2_OnlyChr_.pdf")
pheatmap(assay(norm_Chr2)[Sel_Chr2,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr2, annotation_colors = colors_Chr2, #annotation_row = gff$gene,
         clustering_callback=callback_chr2)
dev.off()

callback_chr4 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr4_Onlychr_.pdf")
pheatmap(assay(norm_Chr4)[Sel_Chr4,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr4, annotation_colors = colors_Chr4, #annotation_row = gff$gene
         clustering_callback=callback_chr4)
dev.off()

callback_chr11 = function(hc, mat){
  n <- length(hc$order)
  sv <- if (n == nrow(mat)) svd(mat)$u[,1] else svd(t(mat))$u[,1]
  dend <- reorder(as.dendrogram(hc), wts = sv[1:n])
  as.hclust(dend)
}
genes_validos11 <- Sel_Chr11[Sel_Chr11 %in% rownames(norm_Chr11)]

pdf(file="heatmap_Chr11_OnlyChr_.pdf")
pheatmap(
  assay(norm_Chr11)[genes_validos11, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  scale = "row",
  annotation_col = annotation_Chr11,
  annotation_colors = colors_Chr11,
  clustering_callback = callback_chr11
)
dev.off()




callback_chr16 = function(hc, mat){
  n <- length(hc$order)
  sv <- if (n == nrow(mat)) svd(mat)$u[,1] else svd(t(mat))$u[,1]
  dend <- reorder(as.dendrogram(hc), wts = sv[1:n])
  as.hclust(dend)
}
genes_validos16 <- Sel_Chr16[Sel_Chr16 %in% rownames(norm_Chr16)]

pdf(file="heatmap_Chr16_OnlyChr_.pdf")
pheatmap(
  assay(norm_Chr16)[genes_validos16, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  scale = "row",
  annotation_col = annotation_Chr16,
  annotation_colors = colors_Chr16,
  clustering_callback = callback_chr16
)
dev.off()



####      Generate PCAs only with data from relevant differentially expessed genes####

PCA_Chr2 <- t(assay(norm_Chr2)[Sel_Chr2,])
PCA_Chr4 <- t(assay(norm_Chr4)[Sel_Chr4,])
PCA_Chr11 <- t(assay(norm_Chr11)[genes_validos11,])
PCA_Chr16 <- t(assay(norm_Chr16)[genes_validos16,])

#do the principal component analyses
PCA_Chr2 <- prcomp(PCA_Chr2)
PCA_Chr4 <- prcomp(PCA_Chr4)
PCA_Chr11 <- prcomp(PCA_Chr11)
PCA_Chr16 <- prcomp(PCA_Chr16)

#get your gene loadings to the PCAs
Chr2_loadings <- as.data.frame(unclass(PCA_Chr2$rotation))
Chr4_loadings <- as.data.frame(unclass(PCA_Chr4$rotation))
Chr11_loadings <- as.data.frame(unclass(PCA_Chr11$rotation))
Chr16_loadings <- as.data.frame(unclass(PCA_Chr16$rotation))

#Pop_loadings %>% arrange(desc(abs(PC1)))
#Get the % expained by each axis
Chr2.eig <- get_eigenvalue(PCA_Chr2)
Chr4.eig <- get_eigenvalue(PCA_Chr4)
Chr11.eig <- get_eigenvalue(PCA_Chr11)
Chr16.eig <- get_eigenvalue(PCA_Chr16)

#generate a matrix to plot
PCA_Chr2 <- as.data.frame(PCA_Chr2$x)
PCA_Chr4 <- as.data.frame(PCA_Chr4$x)
PCA_Chr11 <- as.data.frame(PCA_Chr11$x)
PCA_Chr16 <- as.data.frame(PCA_Chr16$x)

PCA_Chr2 <- cbind(PCA_Chr2, coldata)
PCA_Chr4 <- cbind(PCA_Chr4, coldata4)
PCA_Chr11 <- cbind(PCA_Chr11, coldata11)
PCA_Chr16 <- cbind(PCA_Chr16, coldata)


pdf("PCA_Chr2_genes_onlychr.pdf")
ggplot(data = PCA_Chr2, aes(x = PC1, y = PC2, color=Chr2))+
  geom_segment(data = Chr2_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr2), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr2, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr2.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr2.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr4_genes_onlychr.pdf")
ggplot(data = PCA_Chr4, aes(x = PC1, y = PC2, color=Chr4))+
  geom_segment(data = Chr4_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr4), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr4, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr4.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr4.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr11_genes_onlychr.pdf")
ggplot(data = PCA_Chr11, aes(x = PC1, y = PC2, color=Chr11))+
  geom_segment(data = Chr11_loadings, aes(x = 0, xend = 20*PC1, y = 0, yend = 20*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(color=Chr11, fill=Chr11), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr11, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr11.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr11.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr16_genes_onlychr.pdf")
ggplot(data = PCA_Chr16, aes(x = PC1, y = PC2, color=Chr16))+
  geom_segment(data = Chr16_loadings, aes(x = 0, xend = 5*PC1, y = 0, yend = 5*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr16), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr16, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr16.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr16.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()





#####
#####
#####
####      Differential expression using only genes outside of inversions######
ddseqPop <- DESeq(ddsPop)
ddseqChr2 <- DESeq(ddsChr2)
ddseqChr4 <- DESeq(ddsChr4)
ddseqChr11 <- DESeq(ddsChr11)
ddseqChr16 <- DESeq(ddsChr16)

gene_no_inv <- gff[-which(gff$Chrom=="Chromosome_2" & gff$start>=1 & gff$end<=19480816),]
gene_no_inv <- gene_no_inv[-which(gene_no_inv$Chrom=="Chromosome_4" & gene_no_inv$start>=1 & gene_no_inv$end<=16999389),]
gene_no_inv <- gene_no_inv[-which(gene_no_inv$Chrom=="Chromosome_11" & gene_no_inv$start>=1 & gene_no_inv$end<=10440662),]
gene_no_inv <- gene_no_inv[-which(gene_no_inv$Chrom=="Chromosome_16" & gene_no_inv$start>=4371327 & gene_no_inv$end<=18293200),]

#correct qvalues for only those genes in a  certain chromosome
#install.packages("FDRestimation")
library(FDRestimation)

pval_pop <- as.data.frame(results(ddseqPop))
pval_pop <- pval_pop[gene_no_inv$gene,]
pval_pop <- pval_pop[!is.na(pval_pop$pvalue), ]
pval_pop_test <- pval_pop$pvalue
FDR_pop <-p.fdr(pvalues=pval_pop_test, just.fdr=T)
pval_pop$FDR <- FDR_pop
Sel_Pop <- rownames(subset(pval_pop, FDR < 0.05 & baseMean > 2))

pval_chr2 <- as.data.frame(results(ddseqChr2))
pval_chr2 <- pval_chr2[gene_no_inv$gene,]
pval_chr2 <- pval_chr2[!is.na(pval_chr2$pvalue), ]
pval_chr2_test <- pval_chr2$pvalue
FDR_Chr2 <-p.fdr(pvalues=pval_chr2_test, just.fdr=T)
pval_chr2$FDR <- FDR_Chr2
Sel_Chr2 <- rownames(subset(pval_chr2, FDR < 0.05 & baseMean > 2))

pval_chr4 <- as.data.frame(results(ddseqChr4))
pval_chr4 <- pval_chr4[gene_no_inv$gene,]
pval_chr4 <- pval_chr4[!is.na(pval_chr4$pvalue), ]
pval_chr4_test <- pval_chr4$pvalue
FDR_chr4 <-p.fdr(pvalues=pval_chr4_test, just.fdr=T)
pval_chr4$FDR <- FDR_chr4
Sel_Chr4 <- rownames(subset(pval_chr4, FDR < 0.05 & baseMean > 2))

pval_chr11 <- as.data.frame(results(ddseqChr11))
pval_chr11 <- pval_chr11[gene_no_inv$gene,]
pval_chr11 <- pval_chr11[!is.na(pval_chr11$pvalue), ]
pval_chr11_test <- pval_chr11$pvalue
FDR_chr11 <-p.fdr(pvalues=pval_chr11_test, just.fdr=T)
pval_chr11$FDR <- FDR_chr11
Sel_Chr11 <- rownames(subset(pval_chr11, FDR < 0.05 & baseMean > 2))

pval_chr16 <- as.data.frame(results(ddseqChr16))
pval_chr16 <- pval_chr16[gene_no_inv$gene,]
pval_chr16 <- pval_chr16[!is.na(pval_chr16$pvalue), ]
pval_chr16_test <- pval_chr16$pvalue
FDR_chr16 <-p.fdr(pvalues=pval_chr16_test, just.fdr=T)
pval_chr16$FDR <- FDR_chr16
Sel_Chr16 <- rownames(subset(pval_chr16, FDR < 0.05 & baseMean > 2))

#get an upset plot
list_all <- list(Pop= Sel_Pop, 
                 Chr2= Sel_Chr2, 
                 Chr4= Sel_Chr4,
                 Chr11=Sel_Chr11,
                 Chr16=Sel_Chr16
)

library(ggVennDiagram)
venn_RNAexp<- ggVennDiagram(list_all, label= "c", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none")

ggsave(venn_RNAexp, file="Venn_RNAexp_Noinv.svg", device = "svg", height=5, width=5)

names<-names(list_all)

list_all <- t(splitstackshape:::charMat(list_all, fill = 0, mode="binary"))
colnames(list_all) <- names 
list_all <- rev(as.data.frame(list_all))
sets <-    colnames(list_all) 

pdf("RNA_upset_Noinv.pdf", width=6, height=5.5)
UpSetR::upset(list_all, sets = sets, mb.ratio = c(0.7, 0.3), 
              # sets.bar.color = c("purple","#b05446","#9dcfdd"), 
              order.by = c("freq","degree"), keep.order=T, line.size = 0.8
              ,sets.x.label = "Number of Diff. Expressed RNA", mainbar.y.label = "Number of Diff. Expressed RNA")
dev.off()
dev.off()

#Remove the genes in common (interactions?)
vector1 <- intersect(Sel_Pop, Sel_Chr2)
vector2 <- intersect(Sel_Pop, Sel_Chr4)
vector3 <- intersect(Sel_Pop, Sel_Chr11)
vector4 <- intersect(Sel_Pop, Sel_Chr16)
vector5 <- intersect(Sel_Chr2, Sel_Chr4)
vector6 <- intersect(Sel_Chr2, Sel_Chr11)
vector7 <- intersect(Sel_Chr2, Sel_Chr16)
vector8 <- intersect(Sel_Chr4, Sel_Chr11)
vector9 <- intersect(Sel_Chr4, Sel_Chr16)
vector10 <- intersect(Sel_Chr11, Sel_Chr16)

interact <- c(vector1, vector2, vector3, vector4,vector5,vector6,vector7,vector8,vector9,vector10)
interact <- unique(interact)

#We keep only those genes we know are specific for the factor we tested
Sel_Pop <- setdiff(Sel_Pop, interact)
Sel_Chr2 <- setdiff(Sel_Chr2, interact)
Sel_Chr4 <- setdiff(Sel_Chr4, interact)
Sel_Chr11 <- setdiff(Sel_Chr11, interact)
Sel_Chr16 <- setdiff(Sel_Chr16, interact)

#Export the gene names that are differentially expressed
write.table(Sel_Pop, file="diff_Pop_genes_Noinv.txt", row.names=F, col.names=F)
write.table(Sel_Chr2, file="diff_chr2_genes_Noinv.txt", row.names=F, col.names=F)
write.table(Sel_Chr4, file="diff_chr4_genes_Noinv.txt", row.names=F, col.names=F)
write.table(Sel_Chr11, file="diff_chr11_genes_Noinv.txt", row.names=F, col.names=F)
write.table(Sel_Chr16, file="diff_chr16_genes_Noinv.txt", row.names=F, col.names=F)

####      Plot heatmaps####

#normalize our read matrix 
norm_Pop <- normTransform(ddseqPop)
norm_Chr2 <- normTransform(ddseqChr2)
norm_Chr4 <- normTransform(ddseqChr4)
norm_Chr11 <- normTransform(ddseqChr11)
norm_Chr16 <- normTransform(ddseqChr16)

colors_Pop = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                  Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue",Chromosome_5="black",
                            Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                            Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                            Chromosome_16="lightblue"))
colors_Chr2 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr2 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_2="lightblue"))
colors_Chr4 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr4 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_4="lightblue"))
colors_Chr11 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr11 = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
                    Chrom = c(Chromosome_11="lightblue"))
colors_Chr16 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr16 = c(AA="#ecd6c0", AB="#834f23"),
                    Chrom = c(Chromosome_16="lightblue"))

pdf(file="heatmap_Pop_Noinv.pdf")
pheatmap(assay(norm_Pop)[Sel_Pop,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Pop, annotation_colors = colors_Pop, #annotation_row = gff$gene,
)
dev.off()


pdf(file="heatmap_Chr2_Noinv.pdf")
pheatmap(assay(norm_Chr2)[Sel_Chr2,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr2, annotation_colors = colors_Chr2, #annotation_row = gff$gene,
)
dev.off()

pdf(file="heatmap_Chr4_Noinv.pdf")
pheatmap(assay(norm_Chr4)[Sel_Chr4,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr4, annotation_colors = colors_Chr4, #annotation_row = gff$gene
)
dev.off()

callback_chr11 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_Chr11_Noinv.pdf")
pheatmap(assay(norm_Chr11)[Sel_Chr11,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr11, annotation_colors = colors_Chr11, #annotation_row = gff$gene,
         clustering_callback=callback_chr11)
dev.off()

pdf(file="heatmap_Chr16_Noinv.pdf")
pheatmap(assay(norm_Chr16)[Sel_Chr16,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr16, annotation_colors = colors_Chr16, #annotation_row = gff$gene,
)
dev.off()



####      Generate PCAs only with data from relevant differentially expessed genes####

PCA_Pop <- t(assay(norm_Pop)[Sel_Pop,])
PCA_Chr2 <- t(assay(norm_Chr2)[Sel_Chr2,])
PCA_Chr4 <- t(assay(norm_Chr4)[Sel_Chr4,])
PCA_Chr11 <- t(assay(norm_Chr11)[Sel_Chr11,])
PCA_Chr16 <- t(assay(norm_Chr16)[Sel_Chr16,])

#do the principal component analyses
PCA_Pop <- prcomp(PCA_Pop)
PCA_Chr2 <- prcomp(PCA_Chr2)
PCA_Chr4 <- prcomp(PCA_Chr4)
PCA_Chr11 <- prcomp(PCA_Chr11)
PCA_Chr16 <- prcomp(PCA_Chr16)

#get your gene loadings to the PCAs
Pop_loadings <- as.data.frame(unclass(PCA_Pop$rotation))
Chr2_loadings <- as.data.frame(unclass(PCA_Chr2$rotation))
Chr4_loadings <- as.data.frame(unclass(PCA_Chr4$rotation))
Chr11_loadings <- as.data.frame(unclass(PCA_Chr11$rotation))
Chr16_loadings <- as.data.frame(unclass(PCA_Chr16$rotation))

#Pop_loadings %>% arrange(desc(abs(PC1)))
#Get the % expained by each axis
Pop.eig <- get_eigenvalue(PCA_Pop)
Chr2.eig <- get_eigenvalue(PCA_Chr2)
Chr4.eig <- get_eigenvalue(PCA_Chr4)
Chr11.eig <- get_eigenvalue(PCA_Chr11)
Chr16.eig <- get_eigenvalue(PCA_Chr16)

#generate a matrix to plot
PCA_Pop <- as.data.frame(PCA_Pop$x)
PCA_Chr2 <- as.data.frame(PCA_Chr2$x)
PCA_Chr4 <- as.data.frame(PCA_Chr4$x)
PCA_Chr11 <- as.data.frame(PCA_Chr11$x)
PCA_Chr16 <- as.data.frame(PCA_Chr16$x)

PCA_Pop <- cbind(PCA_Pop, coldata)
PCA_Chr2 <- cbind(PCA_Chr2, coldata)
PCA_Chr4 <- cbind(PCA_Chr4, coldata4)
PCA_Chr11 <- cbind(PCA_Chr11, coldata11)
PCA_Chr16 <- cbind(PCA_Chr16, coldata)




pdf("PCA_Pop_genes_Noinv.pdf")
ggplot(data = PCA_Pop, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Pop_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr2), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Pop, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Pop.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Pop.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr2_genes_Noinv.pdf")
ggplot(data = PCA_Chr2, aes(x = PC1, y = PC2, color=Chr2))+
  geom_segment(data = Chr2_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr2), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr2, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr2.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr2.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr4_genes_Noinv.pdf")
ggplot(data = PCA_Chr4, aes(x = PC1, y = PC2, color=Chr4))+
  geom_segment(data = Chr4_loadings, aes(x = 0, xend = 10*PC1, y = 0, yend = 10*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr4), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr4, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr4.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr4.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr11_genes_Noinv.pdf")
ggplot(data = PCA_Chr11, aes(x = PC1, y = PC2, color=Chr11))+
  geom_segment(data = Chr11_loadings, aes(x = 0, xend = 20*PC1, y = 0, yend = 20*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(color=Chr11, fill=Chr11), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr11, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr11.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr11.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_Chr16_genes_Noinv.pdf")
ggplot(data = PCA_Chr16, aes(x = PC1, y = PC2, color=Chr16))+
  geom_segment(data = Chr16_loadings, aes(x = 0, xend = 5*PC1, y = 0, yend = 5*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Chr16), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr16, size = 4, pch= 19,stroke = 1.5)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  scale_fill_manual(values=c("#ecd6c0", "#834f23",  "#45322E"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr16.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr16.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()


#####
#####
#####
####      Differential expression among populations within each abundant karyotype######

dds2AA <- DESeqDataSetFromMatrix(countData = cts[,-c(1,2,11,13)],
                                  colData = coldata[-c(1,2,11,13),],
                                  design = ~ Pop) #without interactions as we do not have all levels interacting in most of the cases
dds4AB <- DESeqDataSetFromMatrix(countData = cts[,c(1,4,5,11,12,15)],
                                  colData = coldata[c(1,4,5,11,12,15),],
                                  design = ~ Pop)
dds11AB <- DESeqDataSetFromMatrix(countData = cts[,c(6,7,9,11,13,14)],
                                  colData = coldata[c(6,7,9,11,13,14),],
                                  design = ~ Pop)
dds16AA <- DESeqDataSetFromMatrix(countData = cts[,-c(1,3,11:15)],
                                   colData = coldata[-c(1,3,11:15),],
                                   design = ~ Pop)


####      Differential expression (populations and inversions) using all genes######


#set the references for each factor

dds2AA<-dds2AA[-genesToRemove,]
dds4AB<-dds4AB[-genesToRemove,]
dds11AB<-dds11AB[-genesToRemove,]
dds16AA<-dds16AA[-genesToRemove,]


#analyze differential expression!
ddseq2AA <- DESeq(dds2AA)
ddseq4AB <- DESeq(dds4AB)
ddseq11AB <- DESeq(dds11AB)
ddseq16AA <- DESeq(dds16AA)


#set p-value as 0.05 to be significantly different and order according to it
res2AA_95 <- results(ddseq2AA, alpha=0.05)
res4AB_95 <- results(ddseq4AB, alpha=0.05)
res11AB_95 <- results(ddseq11AB, alpha=0.05)
res16AA_95 <- results(ddseq16AA, alpha=0.05)


res2AA_95 <- res2AA_95[order(res2AA_95$pvalue),]
res4AB_95 <- res4AB_95[order(res4AB_95$pvalue),]
res11AB_95 <- res11AB_95[order(res11AB_95$pvalue),]
res16AA_95 <- res16AA_95[order(res16AA_95$pvalue),]

#get summaries and export tables
summary(res2AA_95)
summary(res4AB_95)
summary(res11AB_95)
summary(res16AA_95)


write.table(res2AA_95, file="DESeq_2AA_95.txt", sep="\t", row.names=TRUE)
write.table(res4AB_95, file="DESeq_4AB_95.txt", sep="\t", row.names=TRUE)
write.table(res11AB_95, file="DESeq_11AB_95.txt", sep="\t", row.names=TRUE)
write.table(res16AA_95, file="DESeq_16AA_95.txt", sep="\t", row.names=TRUE)

####      Identify genes expressed differently and represent them ####

#p adjust value below 0.05 and with a mean number of reads of 2
Sig_2AA_95 <- subset(res2AA_95, padj < 0.05 & baseMean > 2) ## NO DIFFERENCES BETWEEN baseMean 2 & 10
Sig_4AB_95 <- subset(res4AB_95, padj < 0.05 & baseMean > 2 )
Sig_11AB_95 <- subset(res11AB_95, padj < 0.05 & baseMean > 2 )
Sig_16AA_95 <- subset(res16AA_95, padj < 0.05 & baseMean > 2 )


summary(Sig_2AA_95)
summary(Sig_4AB_95)
summary(Sig_11AB_95)
summary(Sig_16AA_95)


Sel_2AA <- rownames(Sig_2AA_95)
Sel_4AB <-rownames(Sig_4AB_95)
Sel_11AB <-rownames(Sig_11AB_95)
Sel_16AA <-rownames(Sig_16AA_95)


#get an upset plot
list_all <- list(Pop= Sel_2AA, 
                 Chr2= Sel_4AB, 
                 Chr4= Sel_11AB,
                 Chr11=Sel_16AA
)

library(ggVennDiagram)
venn_RNAexp<- ggVennDiagram(list_all, label= "c", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none")

ggsave(venn_RNAexp, file="Venn_RNAexp_karyotype.svg", device = "svg", height=5, width=5)

names<-names(list_all)

list_all <- t(splitstackshape:::charMat(list_all, fill = 0, mode="binary"))
colnames(list_all) <- names 
list_all <- rev(as.data.frame(list_all))
sets <-    colnames(list_all) 

pdf("RNA_upset_karyotype.pdf", width=6, height=5.5)
UpSetR::upset(list_all, sets = sets, mb.ratio = c(0.7, 0.3), 
              # sets.bar.color = c("purple","#b05446","#9dcfdd"), 
              order.by = c("freq","degree"), keep.order=T, line.size = 0.8
              ,sets.x.label = "Number of Diff. Expressed RNA", mainbar.y.label = "Number of Diff. Expressed RNA")
dev.off()



#Remove the genes in common (interactions?)

vector1 <- intersect(Sel_2AA, Sel_4AB)
vector2 <- intersect(Sel_2AA, Sel_11AB)
vector3 <- intersect(Sel_2AA, Sel_16AA)
vector4 <- intersect(Sel_4AB, Sel_11AB)
vector5 <- intersect(Sel_4AB, Sel_16AA)
vector6 <- intersect(Sel_11AB, Sel_16AA)


interact <- c(vector1, vector2, vector3, vector4,vector5,vector6)
interact <- unique(interact)

#We keep only those genes we know are specific for the factor we tested
Sel_2AA <- setdiff(Sel_2AA, interact)
Sel_4AB <- setdiff(Sel_4AB, interact)
Sel_11AB <- setdiff(Sel_11AB, interact)
Sel_16AA <- setdiff(Sel_16AA, interact)

#Export the gene names that are differentially expressed
write.table(Sel_2AA, file="diff_2AA_genes.txt", row.names=F, col.names=F)
write.table(Sel_4AB, file="diff_4AB_genes.txt", row.names=F, col.names=F)
write.table(Sel_11AB, file="diff_11AB_genes.txt", row.names=F, col.names=F)
write.table(Sel_16AA, file="diff_16AA_genes.txt", row.names=F, col.names=F)

####      Plot heatmaps####

#normalize our read matrix 
norm_2AA <- normTransform(ddseq2AA)
norm_4AB <- normTransform(ddseq4AB)
norm_11AB <- normTransform(ddseq11AB)
norm_16AA <- normTransform(ddseq16AA)


nomral_counts<-assay(norm_2AA)
write.table(nomral_counts, file="normal_counts_2AA.txt", sep="\t", row.names=TRUE)

#generate annotation files for the pheatmap and give them proper rownames
coldata2 <- coldata[-c(1,2,11,13),]
coldata4 <- coldata[c(1,4,5,11,12,15),]
coldata11 <- coldata[c(6,7,9,11,13,14),]
coldata16 <- coldata[-c(1,3,11:15),]

annotation_Chr2 <- data.frame(Chr2=coldata$Chr2, Pop=coldata$Pop)
annotation_Chr4 <- data.frame(Chr4=coldata4$Chr4, Pop=coldata4$Pop)
annotation_Chr11 <- data.frame(Chr11=coldata11$Chr11, Pop=coldata11$Pop)
annotation_Chr16 <- data.frame(Chr16=coldata$Chr16, Pop=coldata$Pop)

rownames(annotation_Chr2) <- rownames(coldata2)
rownames(annotation_Chr4) <- rownames(coldata4) 
rownames(annotation_Chr11) <- rownames(coldata11) 
rownames(annotation_Chr16) <- rownames(coldata16)

#generate color file for our annotation
colors_Chr2 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr2 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                             Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                             Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                             Chromosome_16="lightblue4"))
colors_Chr4 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                   Chr4 = c(AA="#ecd6c0", AB="#834f23"),
                   Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue",Chromosome_5="black",
                             Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                             Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                             Chromosome_16="lightblue4"))
colors_Chr11 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr11 = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
                    Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                              Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                              Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                              Chromosome_16="lightblue4"))
colors_Chr16 = list(Pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
                    Chr16 = c(AA="#ecd6c0", AB="#834f23"),
                    Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
                              Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
                              Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
                              Chromosome_16="lightblue"))


#PLOT. JUST DO IT!!!!

callback_chr2 = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_2AA.pdf", width=7, height =5)
pheatmap(assay(norm_2AA)[Sel_Chr2,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr2, annotation_colors = colors_Chr2, #annotation_row = gff$gene,
         clustering_callback=callback_chr2)
dev.off()


callback_chr4 = function(hc, mat){
  sv = svd(t(mat))$v[,14]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_4AB.pdf", width=7, height =5)
pheatmap(assay(norm_4AB)[Sel_Chr4,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr4, annotation_colors = colors_Chr4, #annotation_row = gff$gene
         )
dev.off()


callback_chr11 = function(hc, mat){
  sv = svd(t(mat))$v[,4]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_11AB.pdf", width=6.5, height =5)
pheatmap(assay(norm_11AB)[Sel_Chr11,], 
         cluster_rows=T, show_rownames=F, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr11, annotation_colors = colors_Chr11, #annotation_row = gff$gene,
         clustering_callback=callback_chr11)
dev.off()


callback_chr16 = function(hc, mat){
  sv = svd(t(mat))$v[,2]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pdf(file="heatmap_16AA.pdf", width=7, height =5)
pheatmap(assay(norm_16AA)[Sel_Chr16,], 
         cluster_rows=T, show_rownames=T, scale="row",
         cluster_cols=T, annotation_col=annotation_Chr16, annotation_colors = colors_Chr16, #annotation_row = gff$gene,
         clustering_callback=callback_chr16)
dev.off()
dev.off()

####      Generate PCAs only with data from relevant differentially expessed genes####

#extract the matrices
PCA_Chr2 <- t(assay(norm_2AA)[Sel_Chr2,])
PCA_Chr4 <- t(assay(norm_4AB)[Sel_Chr4,])
PCA_Chr11 <- t(assay(norm_11AB)[Sel_Chr11,])
PCA_Chr16 <- t(assay(norm_16AA)[Sel_Chr16,])

#do the principal component analyses
PCA_Chr2 <- prcomp(PCA_Chr2)
PCA_Chr4 <- prcomp(PCA_Chr4)
PCA_Chr11 <- prcomp(PCA_Chr11)
PCA_Chr16 <- prcomp(PCA_Chr16)

#get your gene loadings to the PCAs
Chr2_loadings <- as.data.frame(unclass(PCA_Chr2$rotation))
Chr4_loadings <- as.data.frame(unclass(PCA_Chr4$rotation))
Chr11_loadings <- as.data.frame(unclass(PCA_Chr11$rotation))
Chr16_loadings <- as.data.frame(unclass(PCA_Chr16$rotation))

#Pop_loadings %>% arrange(desc(abs(PC1)))
#Get the % expained by each axis
Chr2.eig <- get_eigenvalue(PCA_Chr2)
Chr4.eig <- get_eigenvalue(PCA_Chr4)
Chr11.eig <- get_eigenvalue(PCA_Chr11)
Chr16.eig <- get_eigenvalue(PCA_Chr16)

#generate a matrix to plot
PCA_Chr2 <- as.data.frame(PCA_Chr2$x)
PCA_Chr4 <- as.data.frame(PCA_Chr4$x)
PCA_Chr11 <- as.data.frame(PCA_Chr11$x)
PCA_Chr16 <- as.data.frame(PCA_Chr16$x)

PCA_Chr2 <- cbind(PCA_Chr2, coldata)
PCA_Chr4 <- cbind(PCA_Chr4, coldata4)
PCA_Chr11 <- cbind(PCA_Chr11, coldata11)
PCA_Chr16 <- cbind(PCA_Chr16, coldata)


pdf("PCA_2AA_genes.pdf")
ggplot(data = PCA_Chr2, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Chr2_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr2, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr2.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr2.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()


pdf("PCA_4AB_genes.pdf")
ggplot(data = PCA_Chr4, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Chr4_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr4, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr4.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr4.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()


pdf("PCA_11AB_genes.pdf")
ggplot(data = PCA_Chr11, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Chr11_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr11, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr11.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr11.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()


pdf("PCA_16AA_genes.pdf")
ggplot(data = PCA_Chr16, aes(x = PC1, y = PC2, color=Pop))+
  geom_segment(data = Chr16_loadings, aes(x = 0, xend = 50*PC1, y = 0, yend = 50*PC2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1, alpha=0.4) +
  stat_ellipse(geom="polygon", level=0.9, aes(fill=Pop), alpha=0.1, lwd=1.5)+
  geom_point(data = PCA_Chr16, size = 8, pch= 19,stroke = 1.5, alpha=0.8)+
  # geom_text(data = Pop_loadings[1,], aes(x= 15*PC1, y = 15*PC2, label = rownames(Pop_loadings)), 
  #           size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_manual(values=c("indianred", "steelblue",  "yellow3"))+
  scale_fill_manual(values=c("indianred", "steelblue",  "yellow3"))+
  theme_classic()+
  labs(x = paste("PC1 (", round(Chr16.eig[1,2], 2), "%)", sep = ""), y = paste("PC2 (", round(Chr16.eig[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()





