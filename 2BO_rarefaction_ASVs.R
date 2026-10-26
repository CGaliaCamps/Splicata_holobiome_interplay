library(vegan)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(gplots)
library(car)
library(pairwiseAdonis)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)
library(devtools)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

codes <- read.table("metadata.txt", header = TRUE, sep = "\t", dec = ".")
micro <- read.table("Abundance_matrix_nozero.txt", sep="\t", header=T)
micro <- micro[,c(1:15)]
microt<-t(micro)

# We make a binary matrix
microb<-micro
microb[microb>0]<-1
mean(rowSums(microb))
mean(colSums(microb))
plot(density(rowSums(microb)))
plot(density(colSums(microb)))

min(colSums(micro))
colSums((micro))[head(order(colSums(micro)))]
barplot(colSums(micro))

coloret<-c()
for (i in 1:dim(codes)[1])
{
  if (codes$Population[i]=="BAR") coloret<-c(coloret,"#2E2E2E")
  if (codes$Population[i]=="BLA") coloret<-c(coloret,"#E67E22")
  if (codes$Population[i]=="VIL") coloret<-c(coloret,"#85C1E9")
}


pdf("rarefaction_microbiome.pdf", width=10, height=6)
par(lwd=1.5)
rarecurve(microt,step=100,xlab="Number of reads",
          ylab="Number of ASVs",label=F, 
          cex.lab=1.2,cex=2,lwd=2,
          col=coloret)

legend ("topleft", bty="n",ncol=1, cex=1.2,pt.cex=2.5,lwd=3,
        col=c("#2E2E2E","#E67E22","#85C1E9"),
        c("Barcelona","Blanes","Vilanova"),
        lty=1,
        text.col=c("black"))
dev.off()




