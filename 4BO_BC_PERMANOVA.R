
library(devtools)
library(vegan)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(gplots)
library(car)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)
library(ggpubr)
library(dplyr)
library(reshape)
library(reshape2)
library(data.table)
library(tidyverse)
library(lme4)
library(rsq)
library(emmeans)
library(stringr)
library(pheatmap)
library(gdata)
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
if (!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = TRUE)
  library(rstudioapi)}


#### Get your table ####
setwd(dirname(getActiveDocumentContext()$path))
codes <- read.table("metadata.txt", header = TRUE)
micro <- read.table("Relativefreq_matrix_nozero.txt", sep="\t", header=T, dec = ".")
micro <- micro[,c(1:15)]
micro_chr4 <- micro[,-14]
micro_chr11 <- micro[,-c(8,10)]


#### Calculate your BC distances ####
BC<-vegdist(t(micro))
BC_chr4<-vegdist(t(micro_chr4))
BC_chr11<-vegdist(t(micro_chr11))

#### Perform your PERMANOVAs ####
perm_pop<- adonis2(BC~codes$Pop, 
                   sqrt.dist=F, by="terms", permutations=10000)
perm_pop

perm_chr2<- adonis2(BC~codes$Chr2, 
                   sqrt.dist=F, by="terms", permutations=10000)
perm_chr2

perm_chr4<- adonis2(BC_chr4~codes[-14,]$Chr4, 
                   sqrt.dist=F, by="terms", permutations=10000)
perm_chr4

perm_chr11<- adonis2(BC_chr11~codes[-c(8,10),]$Chr11, 
                   sqrt.dist=F, by="terms", permutations=10000)
perm_chr11

perm_chr16<- adonis2(BC~codes$Chr16, 
                   sqrt.dist=F, by="terms", permutations=10000)
perm_chr16



permdispPop<- betadisper(BC, codes$Pop, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdispPop)
boxplot(permdispPop)

permdispChr2<- betadisper(BC, codes$Chr2, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdispChr2)
boxplot(permdispChr2)

permdispChr4<- betadisper(BC_chr4, codes[-14,]$Chr4, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdispChr4)
boxplot(permdispChr4)

permdispChr11<- betadisper(BC_chr11, codes[-c(8,10),]$Chr11, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdispChr11)
boxplot(permdispChr11)

permdispChr16<- betadisper(BC, codes$Chr16, type="centroid", bias.adjust = F, sqrt.dist = F, add = T)
plot(permdispChr16)
boxplot(permdispChr16)


permutest(betadisper(BC, codes$Pop, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC,codes$Chr2, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC_chr4,codes[-14,]$Chr4, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC_chr11,codes[-c(8,10),]$Chr11, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))
permutest(betadisper(BC,codes$Chr16, type="centroid", bias.adjust = F, sqrt.dist = F, add = T))

# #### Export relevant tables ####
# write.table(perm_all, "permanova_all.txt", sep ="\t", dec=".")
# write.table(pair_BC_all, "pairwise_all.txt", sep ="\t", dec=".")
