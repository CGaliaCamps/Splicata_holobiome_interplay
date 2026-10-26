####Load your libraries and set wd####
library(clusterSim)
library(vegan)
library(pheatmap)
library(dplyr)
library(tidyr)

setwd("E:/Splicata/16S-RNA-DNA/4Mic_vs_RNA")

#### Load your files and filter them####

### METADATA ###
metadata <- read.table(file="metadata.txt", header=T)

### RNA ###
#RNA count file 
RNA_all <- as.data.frame(read.table("Splicata_counts_15.txt", header=T)) 
colnames(RNA_all) <- sub("_15_sorted.bam", "", colnames(RNA_all))
RNA_all <- RNA_all[rowSums(RNA_all)>0,] 
RNA_all <- as.data.frame(RNA_all)

#diff expressed gene lists
genes_pop <- as.data.frame(read.table("diff_pop_genes.txt", header=F)) 
genes_chr2 <- as.data.frame(read.table("diff_chr2_genes.txt", header=F))
genes_chr4 <- as.data.frame(read.table("diff_chr4_genes.txt", header=F))
genes_chr11 <- as.data.frame(read.table("diff_chr11_genes.txt", header=F))
genes_chr16 <- as.data.frame(read.table("diff_chr16_genes.txt", header=F))

#generate count files for diff expressed genes
RNA_pop <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_pop$V1,])
RNA_chr2 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr2$V1,])
RNA_chr4 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr4$V1,])
RNA_chr11 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr11$V1,])
RNA_chr16 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr16$V1,])
  
### MICROBIOME ###
#Microbiome abundance matrix file
Micro_all <- as.data.frame(read.table("Abundance_matrix_nozero.txt", header=T))
Micro_all <- Micro_all[,1:15]

#indval lists
indval_bla <- as.data.frame(read.table("IndvalBLA75_nodups.txt", header=T))
indval_bar <- as.data.frame(read.table("IndvalBAR75_nodups.txt", header=T))
indval_vil <- as.data.frame(read.table("IndvalVIL75_nodups.txt", header=T))
indval_pop <- rbind(indval_bla, indval_bar, indval_vil) 

indval_2AA <- as.data.frame(read.table("Indval2AA75_nodups.txt", header=T))
indval_2AB <- as.data.frame(read.table("Indval2AB75_nodups.txt", header=T))
indval_chr2 <- rbind(indval_2AA, indval_2AB) 

indval_4AA <- as.data.frame(read.table("Indval4AA75_nodups.txt", header=T))
indval_4AB <- as.data.frame(read.table("Indval4AB75_nodups.txt", header=T))
indval_chr4 <- rbind(indval_4AA, indval_4AB) 

indval_11AB <- as.data.frame(read.table("Indval11AB75_nodups.txt", header=T))
indval_11BB <- as.data.frame(read.table("Indval11BB75_nodups.txt", header=T))
indval_chr11 <- rbind(indval_11AB, indval_11BB) 

indval_16AA <- as.data.frame(read.table("Indval16AA75_nodups.txt", header=T))
indval_16AB <- as.data.frame(read.table("Indval16AB75_nodups.txt", header=T))
indval_chr16 <- rbind(indval_16AA, indval_16AB) 

#generate abundance matrix files for indval
Micro_pop <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_pop$ASV,])
Micro_chr2 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr2$ASV,])
Micro_chr4 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr4$ASV,])
Micro_chr11 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr11$ASV,])
Micro_chr16 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr16$ASV,])

#generate annotation files for microbiome (necessary in next steps)
indval_bla$pop <- "BLA"
indval_bar$pop <- "BAR"
indval_vil$pop <- "VIL"
indval_2AA$kar <- "AA"
indval_2AB$kar <- "AB"
indval_4AA$kar <- "AA"
indval_4AB$kar <- "AB"
indval_11AB$kar <- "AB"
indval_11BB$kar <- "BB"
indval_16AA$kar <- "AA"
indval_16AB$kar <- "AB"

annotation_pop <- rbind(indval_bla, indval_bar, indval_vil)
annotation_Chr2 <- rbind(indval_2AA,indval_2AB)
annotation_Chr4 <- rbind(indval_4AA,indval_4AB)
annotation_Chr11 <- rbind(indval_11AB,indval_11BB)
annotation_Chr16 <- rbind(indval_16AA,indval_16AB)

rownames(annotation_pop) <- annotation_pop$ASV
rownames(annotation_Chr2) <- annotation_Chr2$ASV 
rownames(annotation_Chr4) <- annotation_Chr4$ASV 
rownames(annotation_Chr11) <- annotation_Chr11$ASV 
rownames(annotation_Chr16) <- annotation_Chr16$ASV

annotation_pop <- subset.data.frame(annotation_pop, select=c(pop))
annotation_Chr2 <- subset.data.frame(annotation_Chr2, select=c(kar))
annotation_Chr4 <- subset.data.frame(annotation_Chr4, select=c(kar))
annotation_Chr11 <- subset.data.frame(annotation_Chr11, select=c(kar))
annotation_Chr16 <- subset.data.frame(annotation_Chr16, select=c(kar))

#generate annotation files for microbiome (necessary in next steps)
gff<-read.table("splicata_round3_lnc_sno_trna_notracks.gtf", sep="\t")
gff <- gff[which(gff$V3=="transcript"),]
gff <- gff[,c(1,9)]

colnames(gff) <- c("Chrom", "gene")

GOTERM <- gff #Next in the script

gff$gene <- gsub("; .*", "", gff$gene)
gff$gene <- gsub("transcript_id ", "", gff$gene)
gff$gene <- gsub("\\.[^.]*$", "", gff$gene)
gff$gene <- gsub("-.*", "", gff$gene)

gff <- gff[grepl("^Chrom*", gff$Chrom) ,]

gff <- gff[!duplicated(gff),]

rownames(gff) <- gff$gene
annotation_RNA <- subset.data.frame(gff, select=c(Chrom))


### REMOVE UNNECESARY FILES ###
rm(genes_pop, genes_chr2, genes_chr4, genes_chr11, genes_chr16, 
   indval_bla, indval_bar, indval_vil, indval_pop,
   indval_2AA, indval_2AB, indval_chr2, indval_4AA, indval_4AB, indval_chr4,
   indval_11AB, indval_11BB, indval_chr11, indval_16AA, indval_16AB, indval_chr16, gff)


#### Do your correlations ####

#Normalise your data
RNA_pop <- as.matrix(data.Normalization(RNA_pop, type = "n12", normalization = "row"))
Micro_pop <- as.matrix(data.Normalization(Micro_pop, type = "n12", normalization = "row"))

RNA_chr2 <- as.matrix(data.Normalization(RNA_chr2, type = "n12", normalization = "row"))
Micro_chr2 <- as.matrix(data.Normalization(Micro_chr2, type = "n12", normalization = "row"))

RNA_chr4 <- as.matrix(data.Normalization(RNA_chr4, type = "n12", normalization = "row"))
Micro_chr4 <- as.matrix(data.Normalization(Micro_chr4, type = "n12", normalization = "row"))

RNA_chr11 <- as.matrix(data.Normalization(RNA_chr11, type = "n12", normalization = "row"))
Micro_chr11 <- as.matrix(data.Normalization(Micro_chr11, type = "n12", normalization = "row"))

RNA_chr16 <- as.matrix(data.Normalization(RNA_chr16, type = "n12", normalization = "row"))
Micro_chr16 <- as.matrix(data.Normalization(Micro_chr16, type = "n12", normalization = "row"))

# RNA_all <- as.matrix(data.Normalization(RNA_all, type = "n12", normalization = "row"))
# Micro_all <- as.matrix(data.Normalization(Micro_all, type = "n12", normalization = "row"))

# Generate empty matrices to store correlation coefficients
correlation_pop <- matrix(nrow = nrow(RNA_pop), ncol = nrow(Micro_pop))
rownames(correlation_pop) <- rownames(RNA_pop)
colnames(correlation_pop) <- rownames(Micro_pop)

correlation_chr2 <- matrix(nrow = nrow(RNA_chr2), ncol = nrow(Micro_chr2))
rownames(correlation_chr2) <- rownames(RNA_chr2)
colnames(correlation_chr2) <- rownames(Micro_chr2)

correlation_chr4 <- matrix(nrow = nrow(RNA_chr4), ncol = nrow(Micro_chr4))
rownames(correlation_chr4) <- rownames(RNA_chr4)
colnames(correlation_chr4) <- rownames(Micro_chr4)

correlation_chr11 <- matrix(nrow = nrow(RNA_chr11), ncol = nrow(Micro_chr11))
rownames(correlation_chr11) <- rownames(RNA_chr11)
colnames(correlation_chr11) <- rownames(Micro_chr11)

correlation_chr16 <- matrix(nrow = nrow(RNA_chr16), ncol = nrow(Micro_chr16))
rownames(correlation_chr16) <- rownames(RNA_chr16)
colnames(correlation_chr16) <- rownames(Micro_chr16)

# correlation_all <- matrix(nrow = nrow(RNA_all), ncol = nrow(Micro_all))
# rownames(correlation_all) <- rownames(RNA_all)
# colnames(correlation_all) <- rownames(Micro_all)

#Calculate pairwise coeffitient of correlation
for (i in 1:nrow(RNA_pop)) {
  for (j in 1:nrow(Micro_pop)) {
    correlation_pop[i, j] <- cor(RNA_pop[i, ], Micro_pop[j, ])
  }
}

for (i in 1:nrow(RNA_chr2)) {
  for (j in 1:nrow(Micro_chr2)) {
    correlation_chr2[i, j] <- cor(RNA_chr2[i, ], Micro_chr2[j, ])
  }
}

for (i in 1:nrow(RNA_chr4)) {
  for (j in 1:nrow(Micro_chr4)) {
    correlation_chr4[i, j] <- cor(RNA_chr4[i, ], Micro_chr4[j, ])
  }
}

for (i in 1:nrow(RNA_chr11)) {
  for (j in 1:nrow(Micro_chr11)) {
    correlation_chr11[i, j] <- cor(RNA_chr11[i, ], Micro_chr11[j, ])
  }
}

for (i in 1:nrow(RNA_chr16)) {
  for (j in 1:nrow(Micro_chr16)) {
    correlation_chr16[i, j] <- cor(RNA_chr16[i, ], Micro_chr16[j, ])
  }
}
# 
# for (i in 1:nrow(RNA_all)) {
#   for (j in 1:nrow(Micro_all)) {
#     correlation_all[i, j] <- cor(RNA_all[i, ], Micro_all[j, ])
#   }
# }

#calculate coefficient of determination
determin_pop <- (correlation_pop)^2
determin_chr2 <- (correlation_chr2)^2
determin_chr4 <- (correlation_chr4)^2
determin_chr11 <- (correlation_chr11)^2
determin_chr16 <- (correlation_chr16)^2

#check for positive or negative correlations

# transform your matrices to df
correlation_pop <- as.data.frame(correlation_pop)
correlation_chr2 <- as.data.frame(correlation_chr2)
correlation_chr4 <- as.data.frame(correlation_chr4)
correlation_chr11 <- as.data.frame(correlation_chr11)
correlation_chr16 <- as.data.frame(correlation_chr16)

# transform positive correlation values to "+" and negative values to "-" 
correlation_pop[correlation_pop >= 0] <- "NA"
correlation_pop[correlation_pop <= 0] <- "-"
correlation_pop[correlation_pop == "NA"] <- "+"

correlation_chr2[correlation_chr2 >= 0] <- "NA"
correlation_chr2[correlation_chr2 <= 0] <- "-"
correlation_chr2[correlation_chr2 == "NA"] <- "+"

correlation_chr4[correlation_chr4 >= 0] <- "NA"
correlation_chr4[correlation_chr4 <= 0] <- "-"
correlation_chr4[correlation_chr4 == "NA"] <- "+"

correlation_chr11[correlation_chr11 >= 0] <- "NA"
correlation_chr11[correlation_chr11 <= 0] <- "-"
correlation_chr11[correlation_chr11 == "NA"] <- "+"

correlation_chr16[correlation_chr16 >= 0] <- "NA"
correlation_chr16[correlation_chr16 <= 0] <- "-"
correlation_chr16[correlation_chr16 == "NA"] <- "+"


####Plot your results heatmaps####

breaks_corr <- seq(0, 1, length.out=16)


annot_pop <- list(
  pop = c(BAR="indianred", BLA="steelblue", VIL="yellow3"),
  Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue",Chromosome_5="black",
            Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
            Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
            Chromosome_16="lightblue"))

annot_col2 <- list(
  kar = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
  Chrom = c(Chromosome_1="black", Chromosome_2="lightblue",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
            Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
            Chromosome_11="lightblue4",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
            Chromosome_16="lightblue4"))

annot_col11 <- list(
  kar = c(AA="#ecd6c0", AB="#834f23", BB="#45322E"),
  Chrom = c(Chromosome_1="black", Chromosome_2="lightblue4",Chromosome_3="black",Chromosome_4="lightblue4",Chromosome_5="black",
            Chromosome_6="black",Chromosome_7="black",Chromosome_8="black",Chromosome_9="black",Chromosome_10="black",
            Chromosome_11="lightblue",Chromosome_12="black",Chromosome_13="black",Chromosome_14="black",Chromosome_15="black",
            Chromosome_16="lightblue4"))


callback_pop = function(hc, mat){
  n <- length(hc$order)
  sv <- if (n == nrow(mat)) svd(mat)$u[,1] else svd(t(mat))$u[,1]
  dend <- reorder(as.dendrogram(hc), wts = sv[1:n])
  as.hclust(dend)
}
pdf("heatmap_micro_RNA_pop.pdf", width=12)
pheatmap(determin_pop[rowSums(determin_pop > 0.5) > 0, colSums(determin_pop > 0.5) > 0], 
         cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T, breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_pop,
         display_numbers=correlation_pop[rowSums(determin_pop > 0.5) > 0, colSums(determin_pop > 0.5) > 0],
         annotation_row=annotation_RNA, border_color = F, annotation_colors = annot_pop, clustering_callback = callback_pop,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()

pdf("heatmap_micro_RNA_chr2.pdf", width=4.2)
pheatmap(determin_chr2[rowSums(determin_chr2 > 0.5) > 0, colSums(determin_chr2 > 0.5) > 0],
         cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T, breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr2,
         display_numbers=correlation_chr2[rowSums(determin_chr2 > 0.5) > 0, colSums(determin_chr2 > 0.5) > 0],
         annotation_row=annotation_RNA, border_color = F, annotation_colors = annot_col2,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()

 pdf("heatmap_micro_RNA_chr4.pdf")
 pheatmap(determin_chr4[rowSums(determin_chr4 > 0.5) > 0, colSums(determin_chr4 > 0.5) > 0],
          cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T,breaks=breaks_corr,
          legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr4,
          display_numbers=correlation_chr4[rowSums(determin_chr4 > 0.5) > 0, colSums(determin_chr4 > 0.5) > 0],
          annotation_row=annotation_RNA, border_color = F,
          color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
 dev.off()

pdf("heatmap_micro_RNA_chr11.pdf", width=5.5)
pheatmap(determin_chr11[rowSums(determin_chr11 > 0.5) > 0, colSums(determin_chr11 > 0.5) > 0],
         cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T, breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr11,
         display_numbers=correlation_chr11[rowSums(determin_chr11 > 0.5) > 0, colSums(determin_chr11 > 0.5) > 0],
         annotation_row=annotation_RNA,  border_color = F, annotation_colors = annot_col11,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()

 pdf("heatmap_micro_RNA_chr16.pdf")
 pheatmap(determin_chr16[rowSums(determin_chr16 > 0.5) > 0, colSums(determin_chr16 > 0.5) > 0],
          cluster_rows=T, show_rownames=F, cluster_cols=T, show_colnames=F,breaks=breaks_corr,
          legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr16,
          display_numbers=correlation_chr16[rowSums(determin_chr16 > 0.5) > 0, colSums(determin_chr16 > 0.5) > 0],
          annotation_row=annotation_RNA, border_color = F,
          color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
 dev.off()

# pheatmap(correlation_all
#          , cluster_rows=T, show_rownames=F, cluster_cols=T, show_colnames=F,
#          legend_breaks = c(1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1), scale="none",
#          color = colorRampPalette(c("#5B2C6F", "#E8DAEF","white","white", "#A3E4D7", "#117864"))(20))

#### Select only those correlations which may be indicative of interaction ####

#Select interactions above a certain predefined value
pos_pop <- which(determin_pop > 0.5, arr.ind=TRUE)

pos_chr2 <- which(determin_chr2 > 0.5, arr.ind=TRUE)

pos_chr4 <- which(determin_chr4 > 0.5, arr.ind=TRUE)

pos_chr11 <- which(determin_chr11 > 0.5, arr.ind=TRUE)

pos_chr16 <- which(determin_chr16 > 0.5, arr.ind=TRUE)

#Generate interesting correlation file
pos_pop <- data.frame(Row=rownames(determin_pop)[pos_pop[, 1]], 
                        Col=colnames(determin_pop)[pos_pop[, 2]], 
                        Val=determin_pop[pos_pop])

pos_chr2 <- data.frame(Row=rownames(determin_chr2)[pos_chr2[, 1]], 
                        Col=colnames(determin_chr2)[pos_chr2[, 2]], 
                        Val=determin_chr2[pos_chr2])

pos_chr4 <- data.frame(Row=rownames(determin_chr4)[pos_chr4[, 1]], 
                        Col=colnames(determin_chr4)[pos_chr4[, 2]], 
                        Val=determin_chr4[pos_chr4])

pos_chr11 <- data.frame(Row=rownames(determin_chr11)[pos_chr11[, 1]], 
                        Col=colnames(determin_chr11)[pos_chr11[, 2]], 
                        Val=determin_chr11[pos_chr11])

pos_chr16 <- data.frame(Row=rownames(determin_chr16)[pos_chr16[, 1]], 
                      Col=colnames(determin_chr16)[pos_chr16[, 2]], 
                      Val=determin_chr16[pos_chr16])


write.table(pos_pop, file="correlations_pop.txt")
write.table(pos_chr2, file="correlations_chr2.txt")
write.table(pos_chr4, file="correlations_chr4.txt")
write.table(pos_chr11, file="correlations_chr11.txt")
write.table(pos_chr16, file="correlations_chr16.txt")

####Check the functions

GOTERM <- GOTERM[grepl("^Chrom*", GOTERM$Chrom) ,]

GOTERM_ONT <- GOTERM[grepl("Ontology_term ", GOTERM$gene),]
GOTERM_ONT$gene <- gsub(".*Ontology_term ", "", GOTERM_ONT$gene)
GOTERM_ONT$gene <- gsub("; .*", "", GOTERM_ONT$gene)

GOTERM_NAME <- GOTERM[grepl("Ontology_term ", GOTERM$gene),]
GOTERM_NAME$gene <- gsub("transcript_id ", "", GOTERM_NAME$gene)
GOTERM_NAME$gene <- gsub("-.*", "", GOTERM_NAME$gene)

GOTERM_NAME$GO <- GOTERM_ONT$gene
GOTERM <- GOTERM_NAME[!duplicated(GOTERM_NAME),]
GOTERM <- separate_rows(GOTERM, GO, sep = ",")

pos_pops <- pos_pop$Row
pos_chrs2 <- pos_chr2$Row
pos_chrs4 <- pos_chr4$Row
pos_chrs11 <- pos_chr11$Row
pos_chrs16 <- pos_chr16$Row

GOs_pop <- GOTERM[GOTERM$gene %in% pos_pops,]
GOs_chr2 <- GOTERM[GOTERM$gene %in% pos_chrs2,]
GOs_chr4 <- GOTERM[GOTERM$gene %in% pos_chrs4,]
GOs_chr11 <- GOTERM[GOTERM$gene %in% pos_chrs11,]
GOs_chr16 <- GOTERM[GOTERM$gene %in% pos_chrs16,]




#####Repeat everything with the onlychr dataset
#diff expressed gene lists
genes_chr2 <- as.data.frame(read.table("diff_chr2_genes_onlychr.txt", header=F))
genes_chr4 <- as.data.frame(read.table("diff_chr4_genes_onlychr.txt", header=F))
genes_chr11 <- as.data.frame(read.table("diff_chr11_genes_onlychr.txt", header=F))
genes_chr16 <- as.data.frame(read.table("diff_chr16_genes_onlychr.txt", header=F))

#generate count files for diff expressed genes
RNA_chr2 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr2$V1,])
RNA_chr4 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr4$V1,])
RNA_chr11 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr11$V1,])
RNA_chr16 <- as.data.frame(RNA_all[rownames(RNA_all) %in% genes_chr16$V1,])

### MICROBIOME ###
#Microbiome abundance matrix file
Micro_all <- as.data.frame(read.table("Abundance_matrix_nozero.txt", header=T))
Micro_all <- Micro_all[,1:15]

#indval lists
indval_2AA <- as.data.frame(read.table("Indval2AA75_nodups.txt", header=T))
indval_2AB <- as.data.frame(read.table("Indval2AB75_nodups.txt", header=T))
indval_chr2 <- rbind(indval_2AA, indval_2AB) 

indval_4AA <- as.data.frame(read.table("Indval4AA75_nodups.txt", header=T))
indval_4AB <- as.data.frame(read.table("Indval4AB75_nodups.txt", header=T))
indval_chr4 <- rbind(indval_4AA, indval_4AB) 

indval_11AB <- as.data.frame(read.table("Indval11AB75_nodups.txt", header=T))
indval_11BB <- as.data.frame(read.table("Indval11BB75_nodups.txt", header=T))
indval_chr11 <- rbind(indval_11AB, indval_11BB) 

indval_16AA <- as.data.frame(read.table("Indval16AA75_nodups.txt", header=T))
indval_16AB <- as.data.frame(read.table("Indval16AB75_nodups.txt", header=T))
indval_chr16 <- rbind(indval_16AA, indval_16AB) 

#generate abundance matrix files for indval
Micro_chr2 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr2$ASV,])
Micro_chr4 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr4$ASV,])
Micro_chr11 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr11$ASV,])
Micro_chr16 <- as.data.frame(Micro_all[rownames(Micro_all) %in% indval_chr16$ASV,])

#generate annotation files for microbiome (necessary in next steps)
indval_2AA$kar <- "AA"
indval_2AB$kar <- "AB"
indval_4AA$kar <- "AA"
indval_4AB$kar <- "AB"
indval_11AB$kar <- "AB"
indval_11BB$kar <- "BB"
indval_16AA$kar <- "AA"
indval_16AB$kar <- "AB"

annotation_Chr2 <- rbind(indval_2AA,indval_2AB)
annotation_Chr4 <- rbind(indval_4AA,indval_4AB)
annotation_Chr11 <- rbind(indval_11AB,indval_11BB)
annotation_Chr16 <- rbind(indval_16AA,indval_16AB)

rownames(annotation_Chr2) <- annotation_Chr2$ASV 
rownames(annotation_Chr4) <- annotation_Chr4$ASV 
rownames(annotation_Chr11) <- annotation_Chr11$ASV 
rownames(annotation_Chr16) <- annotation_Chr16$ASV

annotation_Chr2 <- subset.data.frame(annotation_Chr2, select=c(kar))
annotation_Chr4 <- subset.data.frame(annotation_Chr4, select=c(kar))
annotation_Chr11 <- subset.data.frame(annotation_Chr11, select=c(kar))
annotation_Chr16 <- subset.data.frame(annotation_Chr16, select=c(kar))

#generate annotation files for microbiome (necessary in next steps)
gff<-read.table("splicata_round3_lnc_sno_trna_notracks.gtf", sep="\t")
gff <- gff[which(gff$V3=="transcript"),]
gff <- gff[,c(1,9)]

colnames(gff) <- c("Chrom", "gene")

GOTERM <- gff #Next in the script

gff$gene <- gsub("; .*", "", gff$gene)
gff$gene <- gsub("transcript_id ", "", gff$gene)
gff$gene <- gsub("\\.[^.]*$", "", gff$gene)
gff$gene <- gsub("-.*", "", gff$gene)

gff <- gff[grepl("^Chrom*", gff$Chrom) ,]

gff <- gff[!duplicated(gff),]

rownames(gff) <- gff$gene
annotation_RNA <- subset.data.frame(gff, select=c(Chrom))


### REMOVE UNNECESARY FILES ###
rm(genes_pop, genes_chr2, genes_chr4, genes_chr11, genes_chr16, 
   indval_bla, indval_bar, indval_vil, indval_pop,
   indval_2AA, indval_2AB, indval_chr2, indval_4AA, indval_4AB, indval_chr4,
   indval_11AB, indval_11BB, indval_chr11, indval_16AA, indval_16AB, indval_chr16, gff)


#### Do your correlations ####

#Normalise your data
RNA_chr2 <- as.matrix(data.Normalization(RNA_chr2, type = "n12", normalization = "row"))
Micro_chr2 <- as.matrix(data.Normalization(Micro_chr2, type = "n12", normalization = "row"))

RNA_chr4 <- as.matrix(data.Normalization(RNA_chr4, type = "n12", normalization = "row"))
Micro_chr4 <- as.matrix(data.Normalization(Micro_chr4, type = "n12", normalization = "row"))

RNA_chr11 <- as.matrix(data.Normalization(RNA_chr11, type = "n12", normalization = "row"))
Micro_chr11 <- as.matrix(data.Normalization(Micro_chr11, type = "n12", normalization = "row"))

RNA_chr16 <- as.matrix(data.Normalization(RNA_chr16, type = "n12", normalization = "row"))
Micro_chr16 <- as.matrix(data.Normalization(Micro_chr16, type = "n12", normalization = "row"))

# RNA_all <- as.matrix(data.Normalization(RNA_all, type = "n12", normalization = "row"))
# Micro_all <- as.matrix(data.Normalization(Micro_all, type = "n12", normalization = "row"))

# Generate empty matrices to store correlation coefficients
correlation_chr2 <- matrix(nrow = nrow(RNA_chr2), ncol = nrow(Micro_chr2))
rownames(correlation_chr2) <- rownames(RNA_chr2)
colnames(correlation_chr2) <- rownames(Micro_chr2)

correlation_chr4 <- matrix(nrow = nrow(RNA_chr4), ncol = nrow(Micro_chr4))
rownames(correlation_chr4) <- rownames(RNA_chr4)
colnames(correlation_chr4) <- rownames(Micro_chr4)

correlation_chr11 <- matrix(nrow = nrow(RNA_chr11), ncol = nrow(Micro_chr11))
rownames(correlation_chr11) <- rownames(RNA_chr11)
colnames(correlation_chr11) <- rownames(Micro_chr11)

correlation_chr16 <- matrix(nrow = nrow(RNA_chr16), ncol = nrow(Micro_chr16))
rownames(correlation_chr16) <- rownames(RNA_chr16)
colnames(correlation_chr16) <- rownames(Micro_chr16)

# correlation_all <- matrix(nrow = nrow(RNA_all), ncol = nrow(Micro_all))
# rownames(correlation_all) <- rownames(RNA_all)
# colnames(correlation_all) <- rownames(Micro_all)

#Calculate pairwise coeffitient of correlation
for (i in 1:nrow(RNA_chr2)) {
  for (j in 1:nrow(Micro_chr2)) {
    correlation_chr2[i, j] <- cor(RNA_chr2[i, ], Micro_chr2[j, ])
  }
}

for (i in 1:nrow(RNA_chr4)) {
  for (j in 1:nrow(Micro_chr4)) {
    correlation_chr4[i, j] <- cor(RNA_chr4[i, ], Micro_chr4[j, ])
  }
}

for (i in 1:nrow(RNA_chr11)) {
  for (j in 1:nrow(Micro_chr11)) {
    correlation_chr11[i, j] <- cor(RNA_chr11[i, ], Micro_chr11[j, ])
  }
}

for (i in 1:nrow(RNA_chr16)) {
  for (j in 1:nrow(Micro_chr16)) {
    correlation_chr16[i, j] <- cor(RNA_chr16[i, ], Micro_chr16[j, ])
  }
}
# 
# for (i in 1:nrow(RNA_all)) {
#   for (j in 1:nrow(Micro_all)) {
#     correlation_all[i, j] <- cor(RNA_all[i, ], Micro_all[j, ])
#   }
# }

#calculate coefficient of determination
determin_chr2 <- (correlation_chr2)^2
determin_chr4 <- (correlation_chr4)^2
determin_chr11 <- (correlation_chr11)^2
determin_chr16 <- (correlation_chr16)^2

#check for positive or negative correlations

# transform your matrices to df
correlation_chr2 <- as.data.frame(correlation_chr2)
correlation_chr4 <- as.data.frame(correlation_chr4)
correlation_chr11 <- as.data.frame(correlation_chr11)
correlation_chr16 <- as.data.frame(correlation_chr16)

# transform positive correlation values to "+" and negative values to "-" 
correlation_chr2[correlation_chr2 >= 0] <- "NA"
correlation_chr2[correlation_chr2 <= 0] <- "-"
correlation_chr2[correlation_chr2 == "NA"] <- "+"

correlation_chr4[correlation_chr4 >= 0] <- "NA"
correlation_chr4[correlation_chr4 <= 0] <- "-"
correlation_chr4[correlation_chr4 == "NA"] <- "+"

correlation_chr11[correlation_chr11 >= 0] <- "NA"
correlation_chr11[correlation_chr11 <= 0] <- "-"
correlation_chr11[correlation_chr11 == "NA"] <- "+"

correlation_chr16[correlation_chr16 >= 0] <- "NA"
correlation_chr16[correlation_chr16 <= 0] <- "-"
correlation_chr16[correlation_chr16 == "NA"] <- "+"


####Plot your results heatmaps####

breaks_corr <- seq(0, 1, length.out=16)


pdf("heatmap_micro_RNA_chr2_onlychr.pdf", width=6)
pheatmap(determin_chr2[rowSums(determin_chr2 > 0.5) > 0, colSums(determin_chr2 > 0.5) > 0],
         cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T, breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr2,
         display_numbers=correlation_chr2[rowSums(determin_chr2 > 0.5) > 0, colSums(determin_chr2 > 0.5) > 0],
         annotation_row=annotation_RNA, border_color = F, annotation_colors = annot_col2,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()

 pdf("heatmap_micro_RNA_chr4_onlychr.pdf")
 pheatmap(determin_chr4[rowSums(determin_chr4 > 0.5) > 0, colSums(determin_chr4 > 0.5) > 0],
          cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T,breaks=breaks_corr,
          legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr4,
          display_numbers=correlation_chr4[rowSums(determin_chr4 > 0.5) > 0, colSums(determin_chr4 > 0.5) > 0],
          annotation_row=annotation_RNA, border_color = F, annotation_colors = annot_col4,
          color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
 dev.off()

pdf("heatmap_micro_RNA_chr11_onlychr.pdf", width=12)
pheatmap(determin_chr11[rowSums(determin_chr11 > 0.5) > 0, colSums(determin_chr11 > 0.5) > 0],
         cluster_rows=T, show_rownames=T, cluster_cols=T, show_colnames=T, breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr11,
         display_numbers=correlation_chr11[rowSums(determin_chr11 > 0.5) > 0, colSums(determin_chr11 > 0.5) > 0],
         annotation_row=annotation_RNA,  border_color = F, annotation_colors = annot_col11,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()

pdf("heatmap_micro_RNA_chr16_onlychr.pdf")
pheatmap(determin_chr16[rowSums(determin_chr16 > 0.5) > 0, colSums(determin_chr16 > 0.5) > 0],
         cluster_rows=T, show_rownames=F, cluster_cols=T, show_colnames=F,breaks=breaks_corr,
         legend_breaks = c(1,0.8,0.6,0.4,0.2,0), scale="none", annotation_col = annotation_Chr16,
         display_numbers=correlation_chr16[rowSums(determin_chr16 > 0.5) > 0, colSums(determin_chr16 > 0.5) > 0],
         annotation_row=annotation_RNA, border_color = F, annotation_colors = annot_col16,
         color = colorRampPalette(c("white","white", "#D2B4DE", "#633974"))(16))
dev.off()


#### Select only those correlations which may be indicative of interaction ####

#Select interactions above a certain predefined value
pos_chr2 <- which(determin_chr2 > 0.5, arr.ind=TRUE)

pos_chr4 <- which(determin_chr4 > 0.5, arr.ind=TRUE)

pos_chr11 <- which(determin_chr11 > 0.5, arr.ind=TRUE)

pos_chr16 <- which(determin_chr16 > 0.5, arr.ind=TRUE)

#Generate interesting correlation file
pos_chr2 <- data.frame(Row=rownames(determin_chr2)[pos_chr2[, 1]], 
                       Col=colnames(determin_chr2)[pos_chr2[, 2]], 
                       Val=determin_chr2[pos_chr2])

pos_chr4 <- data.frame(Row=rownames(determin_chr4)[pos_chr4[, 1]], 
                       Col=colnames(determin_chr4)[pos_chr4[, 2]], 
                       Val=determin_chr4[pos_chr4])

pos_chr11 <- data.frame(Row=rownames(determin_chr11)[pos_chr11[, 1]], 
                        Col=colnames(determin_chr11)[pos_chr11[, 2]], 
                        Val=determin_chr11[pos_chr11])

pos_chr16 <- data.frame(Row=rownames(determin_chr16)[pos_chr16[, 1]], 
                        Col=colnames(determin_chr16)[pos_chr16[, 2]], 
                        Val=determin_chr16[pos_chr16])


write.table(pos_chr2, file="correlations_chr2_onlychr.txt")
write.table(pos_chr4, file="correlations_chr4_onlychr.txt")
write.table(pos_chr11, file="correlations_chr11_onlychr.txt")
write.table(pos_chr16, file="correlations_chr16_onlychr.txt")

####Check the functions

GOTERM <- GOTERM[grepl("^Chrom*", GOTERM$Chrom) ,]

GOTERM_ONT <- GOTERM[grepl("Ontology_term ", GOTERM$gene),]
GOTERM_ONT$gene <- gsub(".*Ontology_term ", "", GOTERM_ONT$gene)
GOTERM_ONT$gene <- gsub("; .*", "", GOTERM_ONT$gene)

GOTERM_NAME <- GOTERM[grepl("Ontology_term ", GOTERM$gene),]
GOTERM_NAME$gene <- gsub("transcript_id ", "", GOTERM_NAME$gene)
GOTERM_NAME$gene <- gsub("-.*", "", GOTERM_NAME$gene)

GOTERM_NAME$GO <- GOTERM_ONT$gene
GOTERM <- GOTERM_NAME[!duplicated(GOTERM_NAME),]
GOTERM <- separate_rows(GOTERM, GO, sep = ",")

pos_chrs2 <- pos_chr2$Row
pos_chrs4 <- pos_chr4$Row
pos_chrs11 <- pos_chr11$Row
pos_chrs16 <- pos_chr16$Row

GOs_chr2 <- GOTERM[GOTERM$gene %in% pos_chrs2,]
GOs_chr4 <- GOTERM[GOTERM$gene %in% pos_chrs4,]
GOs_chr11 <- GOTERM[GOTERM$gene %in% pos_chrs11,]
GOs_chr16 <- GOTERM[GOTERM$gene %in% pos_chrs16,]

