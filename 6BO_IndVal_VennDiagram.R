## Comencem amb les llibreries i el directori
if (!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = TRUE)
  library(rstudioapi)}
if (!require("indicspecies")){
  install.packages("indicspecies", dependencies = TRUE)
  library(indicspecies)}
if (!require("readxl")){
  install.packages("readxl", dependencies = TRUE)
  library(readxl)} 
library(VennDiagram)

##Read abundance table
setwd("E:/Splicata/16S-RNA-DNA/3Microbioma")
sample <-read.table("Abundance_matrix_nozero.txt",stringsAsFactors = F, header=T, sep="\t")
# Remove non-numeric columns
sample <- sample[,1:15]

# Relative Frequencies
sample<-as.matrix(sample)
sample<-prop.table(sample,2)
sample <- as.data.frame(sample)

# Metadata time
codes <- read.table("E:/Splicata/16S-RNA-DNA/metadata.txt",header=T,stringsAsFactors = F)
rownames(codes)<-codes$Sample
codes <- codes[order(codes$Sample),]

#transpose
sample<-t(sample)

group_pop<-factor(codes$Pop,levels=unique(codes$Pop))
group_chr2<-factor(codes$Chr2,levels=unique(codes$Chr2))
group_chr4<-factor(codes$Chr4[-14],levels=unique(codes$Chr4[-14]))
group_chr11<-factor(codes$Chr11[-c(8,10)],levels=unique(codes$Chr11[-c(8,10)]))
group_chr16<-factor(codes$Chr16,levels=unique(codes$Chr16))

indval_pop <- multipatt(sample,group_pop,control=how(nperm=999),max.order=1, duleg=TRUE)
indval_chr2 <- multipatt(sample,group_chr2,control=how(nperm=999),max.order=1, duleg=TRUE)
indval_chr4 <- multipatt(sample[-15,],group_chr4,control=how(nperm=999),max.order=1, duleg=TRUE)
indval_chr11 <- multipatt(sample[-c(9,11),],group_chr11,control=how(nperm=999),max.order=1, duleg=TRUE)
indval_chr16<- multipatt(sample,group_chr16,control=how(nperm=999),max.order=1, duleg=TRUE)



 #Export tables
  sink("indval_pop.txt")
  cat(summary(indval_pop,indvalcomp = T))
  sink()

  sink("indval_chr2.txt")
  cat(summary(indval_chr2,indvalcomp = T))
  sink()
  
  sink("indval_chr4.txt")
  cat(summary(indval_chr4,indvalcomp = T))
  sink()
  
  sink("indval_chr11.txt")
  cat(summary(indval_chr11,indvalcomp = T))
  sink()
  
  sink("indval_chr16.txt")
  cat(summary(indval_chr16,indvalcomp = T))
  sink()




#add taxononomy to the matrix
taxon <- read.table("Abundance_matrix_nozero.txt",stringsAsFactors = F, header=T, sep="\t")

# taxo <- read.csv("taxonomia_final.csv",header=T,sep=",", stringsAsFactors = F)
# #taxo <- gsub("indet","Unknown", taxo)
# rownames(taxon) <- taxon$ASV.ID
# taxo<- taxo[,c(3:10)]
# sample$ID <- taxon$data.ID

#Load Indval data and keep only thos ASV with A & B above 0.75. (75%  reads in 75% samples)

IndvalBAR <- read.table("indval_BAR.txt", header=T, sep="\t")
IndvalBLA <- read.table("indval_BLA.txt", header=T, sep="\t")
IndvalVIL <- read.table("indval_VIL.txt", header=T, sep="\t")
Indval2AA <- read.table("indval_Chr2_AA.txt", header=T, sep="\t")
Indval2AB <- read.table("indval_Chr2_AB.txt", header=T, sep="\t")
Indval4AA <- read.table("indval_Chr4_AA.txt", header=T, sep="\t")
Indval4AB <- read.table("indval_Chr4_AB.txt", header=T, sep="\t")
Indval11AB <- read.table("indval_Chr11_AB.txt", header=T, sep="\t")
Indval11BB <- read.table("indval_Chr11_BB.txt", header=T, sep="\t")
Indval16AA <- read.table("indval_Chr16_AA.txt", header=T, sep="\t")
Indval16AB <- read.table("indval_Chr16_AB.txt", header=T, sep="\t")

IndvalBAR8 <-IndvalBAR[which(IndvalBAR$A>=0.75 & IndvalBAR$B>=0.75),] 
IndvalBLA8 <-IndvalBLA[which(IndvalBLA$A>=0.75 & IndvalBLA$B>=0.75),] 
IndvalVIL8 <-IndvalVIL[which(IndvalVIL$A>=0.75 & IndvalVIL$B>=0.75),] 
Indval2AA8 <-Indval2AA[which(Indval2AA$A>=0.75 & Indval2AA$B>=0.75),] 
Indval2AB8 <-Indval2AB[which(Indval2AB$A>=0.75 & Indval2AB$B>=0.75),] 
Indval4AA8 <-Indval4AA[which(Indval4AA$A>=0.75 & Indval4AA$B>=0.75),] 
Indval4AB8 <-Indval4AB[which(Indval4AB$A>=0.75 & Indval4AB$B>=0.75),] 
Indval11AB8 <-Indval11AB[which(Indval11AB$A>=0.75 & Indval11AB$B>=0.75),] 
Indval11BB8 <-Indval11BB[which(Indval11BB$A>=0.75 & Indval11BB$B>=0.75),] 
Indval16AA8 <-Indval16AA[which(Indval16AA$A>=0.75 & Indval16AA$B>=0.75),] 
Indval16AB8 <-Indval16AB[which(Indval16AB$A>=0.75 & Indval16AB$B>=0.75),] 

#Plot your upset plots
library(UpSetR)

listBAR8 <- IndvalBAR8$ASV
listBLA8 <- IndvalBLA8$ASV
listVIL8 <- IndvalVIL8$ASV
list2AA8 <- Indval2AA8$ASV
list2AB8 <- Indval2AB8$ASV
list4AA8 <- Indval4AA8$ASV
list4AB8 <- Indval4AB8$ASV
list11AB8 <- Indval11AB8$ASV
list11BB8 <- Indval11BB8$ASV
list16AA8 <- Indval16AA8$ASV
list16AB8 <- Indval16AB8$ASV

list_all <- list(BLA= listBLA8, 
                   BAR= listBAR8, 
                   VIL=listVIL8,
                   Chr2_AA=list2AA8,
                   Chr2_AB=list2AB8,
                   Chr4_AA=list4AA8,
                   Chr4_AB=list4AB8,
                   Chr11_AA=list11AB8,
                   Chr11_BB=list11BB8,
                   Chr16_AA=list16AA8,
                   Chr16_AB=list16AB8)

names<-names(list_all)

list_all <- t(splitstackshape:::charMat(list_all, fill = 0, mode="binary"))
colnames(list_all) <- names 
list_all <- rev(as.data.frame(list_all))
sets <-    colnames(list_all) 

pdf("all_correlated.pdf", width=5.5, height=5.5)
UpSetR::upset(list_all, sets = sets, mb.ratio = c(0.75, 0.25), 
              #sets.bar.color = c("purple","#b05446","#9dcfdd"), 
              order.by = c("freq","degree"), keep.order=T, line.size = 0.8
              ,sets.x.label = "Number of ASVs", mainbar.y.label = "Shared ASV")
dev.off()

library(ggplot2)
library(ggVennDiagram)


Venn_plot = list(Pop=c(listBLA8, listBAR8, listVIL8), 
         Chr2=c(list2AA8, list2AB8), 
         Chr4=c(list4AA8, list4AB8),
         Chr11=c(list11AB8, list11BB8),
         Chr16=c(list16AA8, list16AB8))  


venn_micro<- ggVennDiagram(Venn_plot, label= "c", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme(legend.position = "none")

ggsave(venn_micro, file="Venn_micro.svg", device = "svg", height=5, width=5)


write.table(list_all, "recount_Indval.txt")

remove<-list_all[rowSums(list_all) >= 2, , drop = FALSE]

IndvalBAR8 <- IndvalBAR8[!IndvalBAR8$ASV %in% rownames(remove),]
IndvalBLA8 <- IndvalBLA8[!IndvalBLA8$ASV %in% rownames(remove),]
IndvalVIL8 <- IndvalVIL8[!IndvalVIL8$ASV %in% rownames(remove),]
Indval2AA8 <- Indval2AA8[!Indval2AA8$ASV %in% rownames(remove),]
Indval2AB8 <- Indval2AB8[!Indval2AB8$ASV %in% rownames(remove),]
Indval4AA8 <- Indval4AA8[!Indval4AA8$ASV %in% rownames(remove),]
Indval4AB8 <- Indval4AB8[!Indval4AB8$ASV %in% rownames(remove),]
Indval11AB8 <- Indval11AB8[!Indval11AB8$ASV %in% rownames(remove),]
Indval11BB8 <- Indval11BB8[!Indval11BB8$ASV %in% rownames(remove),]
Indval16AA8 <- Indval16AA8[!Indval16AA8$ASV %in% rownames(remove),]
Indval16AB8 <- Indval16AB8[!Indval16AB8$ASV %in% rownames(remove),]

write.table(IndvalBAR8, "IndvalBAR75_nodups.txt")
write.table(IndvalBLA8, "IndvalBLA75_nodups.txt")
write.table(IndvalVIL8, "IndvalVIL75_nodups.txt")
write.table(Indval2AA8, "Indval2AA75_nodups.txt")
write.table(Indval2AB8, "Indval2AB75_nodups.txt")
write.table(Indval4AA8, "Indval4AA75_nodups.txt")
write.table(Indval4AB8, "Indval4AB75_nodups.txt")
write.table(Indval11AB8, "Indval11AB75_nodups.txt")
write.table(Indval11BB8, "Indval11BB75_nodups.txt")
write.table(Indval16AA8, "Indval16AA75_nodups.txt")
write.table(Indval16AB8, "Indval16AB75_nodups.txt")

donut_pop <- data.frame(Population=c("BAR", "BLA", "VIL"), count=c(nrow(IndvalBAR8), nrow(IndvalBLA8), nrow(IndvalVIL8)))
donut_chr2 <- data.frame(Karyotype=c("AA", "AB"), count=c(nrow(Indval2AA8), nrow(Indval2AB8)))
donut_chr4 <- data.frame(Karyotype=c("AA", "AB"), count=c(nrow(Indval4AA8), nrow(Indval4AB8)))
donut_chr11 <- data.frame(Karyotype=c("AB", "BB"), count=c(nrow(Indval11AB8), nrow(Indval11BB8)))
donut_chr16 <- data.frame(Karyotype=c("AA", "AB"), count=c(nrow(Indval16AA8), nrow(Indval16AB8)))


donut_pop$fraction = donut_pop$count / sum(donut_pop$count)
donut_pop$ymax = cumsum(donut_pop$fraction)
donut_pop$ymin = c(0, head(donut_pop$ymax, n=-1))

donut_chr2$fraction = donut_chr2$count / sum(donut_chr2$count)
donut_chr2$ymax = cumsum(donut_chr2$fraction)
donut_chr2$ymin = c(0, head(donut_chr2$ymax, n=-1))

donut_chr4$fraction = donut_chr4$count / sum(donut_chr4$count)
donut_chr4$ymax = cumsum(donut_chr4$fraction)
donut_chr4$ymin = c(0, head(donut_chr4$ymax, n=-1))

donut_chr11$fraction = donut_chr11$count / sum(donut_chr11$count)
donut_chr11$ymax = cumsum(donut_chr11$fraction)
donut_chr11$ymin = c(0, head(donut_chr11$ymax, n=-1))

donut_chr16$fraction = donut_chr16$count / sum(donut_chr16$count)
donut_chr16$ymax = cumsum(donut_chr16$fraction)
donut_chr16$ymin = c(0, head(donut_chr16$ymax, n=-1))


pdf("donut_pop.pdf", width=2, height=2)
ggplot(donut_pop, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Population, color="black")) +
  geom_rect(lwd=0.5) +
  geom_label( x=3.5, aes(y=(donut_pop$ymax + donut_pop$ymin) / 2, label= donut_pop$count), size=4, label.size=0, label.padding=unit(0, "lines")) +
  scale_fill_manual(values=c("indianred","steelblue","yellow3")) +
  scale_color_manual(values=c("black","black","black"))+
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

pdf("donut_chr2.pdf", width=2, height=2)
ggplot(donut_chr2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Karyotype, color="black")) +
  geom_rect(lwd=0.5) +
  geom_label( x=3.5, aes(y=(donut_chr2$ymax + donut_chr2$ymin) / 2, label= donut_chr2$count), size=4, label.size=0, label.padding=unit(0, "lines")) +
  scale_fill_manual(values=c("#ecd6c0","#834f23")) +
  scale_color_manual(values=c("black","black"))+
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

pdf("donut_chr4.pdf", width=2, height=2)
ggplot(donut_chr4, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Karyotype, color="black")) +
  geom_rect(lwd=0.5) +
  geom_label( x=3.5, aes(y=(donut_chr4$ymax + donut_chr4$ymin) / 2, label= donut_chr4$count), size=4, label.size=0, label.padding=unit(0, "lines")) +
  scale_fill_manual(values=c("#ecd6c0","#834f23")) +
  scale_color_manual(values=c("black","black"))+
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

pdf("donut_chr11.pdf", width=2, height=2)
ggplot(donut_chr11, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Karyotype, color="black")) +
  geom_rect(lwd=0.5) +
  geom_label( x=3.5, aes(y=(donut_chr11$ymax + donut_chr11$ymin) / 2, label= donut_chr11$count), size=4, label.size=0, label.padding=unit(0, "lines")) +
  scale_fill_manual(values=c("#834f23","#45322E")) +
  scale_color_manual(values=c("black","black"))+
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

pdf("donut_chr16.pdf", width=2, height=2)
ggplot(donut_chr16, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Karyotype, color="black")) +
  geom_rect(lwd=0.5) +
  geom_label( x=3.5, aes(y=(donut_chr16$ymax + donut_chr16$ymin) / 2, label= donut_chr16$count), size=4, label.size=0, label.padding=unit(0, "lines")) +
  scale_fill_manual(values=c("#ecd6c0","#834f23")) +
  scale_color_manual(values=c("black","black"))+
  coord_polar(theta="y") +
  xlim(c(1, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()
