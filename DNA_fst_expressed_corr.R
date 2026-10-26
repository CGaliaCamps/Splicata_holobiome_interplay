library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")
library("dplyr")
library("scales")
library("ggplot2")
library("plyr")

setwd("E:/Splicata/16S-RNA-DNA/5DNA_vs_RNA")

####Prepare your gff####
gff <- read.table(file="splicata_round3_lnc_sno_trna_notracks.gtf", sep = "\t")
gff <- gff[which(gff$V3=="transcript"),]
gff <- gff[,c(1,4,5,9)]
colnames(gff) <- c("Chrom", "Start", "Finish", "gene")

gff$gene <- gsub("; .*", "", gff$gene)
gff$gene <- gsub("transcript_id ", "", gff$gene)
gff$gene <- gsub("\\.[^.]*$", "", gff$gene)
gff$gene <- gsub("-.*", "", gff$gene)
gff <- gff[grepl("^Chrom*", gff$Chrom) ,]
gff <- gff[!duplicated(gff),]

####Prepare your correlation files####
Cpop <- read.table(file="correlations_pop.txt", sep = " ", header=T)
Cchr2 <- read.table(file="correlations_chr2.txt", sep = " ", header=T)
Cchr2o <- read.table(file="correlations_chr2_onlychr.txt", sep = " ", header=T)
Cchr4 <- read.table(file="correlations_chr4.txt", sep = " ", header=T)
Cchr4o <- read.table(file="correlations_chr4_onlychr.txt", sep = " ", header=T)
Cchr11 <- read.table(file="correlations_chr11.txt", sep = " ", header=T)
Cchr11o <- read.table(file="correlations_chr11_onlychr.txt", sep = " ", header=T)
Cchr16 <- read.table(file="correlations_chr16.txt", sep = " ", header=T)
Cchr16o <- read.table(file="correlations_chr16_onlychr.txt", sep = " ", header=T)

Corr <- rbind(Cpop,Cchr2,Cchr2o,Cchr4,Cchr4o,Cchr11,Cchr11o,Cchr16,Cchr16o)
Corr <- as.data.frame(Corr$Row)
Corr <- unique(Corr)

####Prepare your differential expression files####
Dpop <- read.table(file="diff_pop_genes.txt", sep = " ")
Dchr2 <- read.table(file="diff_chr2_genes.txt", sep = " ")
Dchr2o <- read.table(file="diff_chr2_genes_onlychr.txt", sep = " ")
Dchr4 <- read.table(file="diff_chr4_genes.txt", sep = " ")
Dchr4o <- read.table(file="diff_chr4_genes_onlychr.txt", sep = " ")
Dchr11 <- read.table(file="diff_chr11_genes.txt", sep = " ")
Dchr11o <- read.table(file="diff_chr11_genes_onlychr.txt", sep = " ")
Dchr16 <- read.table(file="diff_chr16_genes.txt", sep = " ")
Dchr16o <- read.table(file="diff_chr16_genes_onlychr.txt", sep = " ")

Diff <- rbind(Dpop,Dchr2,Dchr2o,Dchr4,Dchr4o,Dchr11,Dchr11o,Dchr16,Cchr16o)
Diff <- unique(Diff)

####Prepare your fst files####

df1pop <- read.table(file="tot_100_globalfst_2500.windowed.weir.fst", 
                  sep = "\t", header=TRUE)

df1chr2 <- read.table(file="chr2_AAvsAB_fst2500.windowed.weir.fst", 
                  sep = "\t", header=TRUE)

df1chr4 <- read.table(file="chr4_AAvsAB_fst2500.windowed.weir.fst", 
                  sep = "\t", header=TRUE)

# df1chr11G <- read.table(file="chr11_AAvsABvsBB_fst2500.windowed.weir.fst", 
#                   sep = "\t", header=TRUE)
df1chr11P <- read.table(file="chr11_BBvsAB_fst2500.windowed.weir.fst", 
                  sep = "\t", header=TRUE)

df1chr16 <- read.table(file="chr16_AAvsAB_fst2500.windowed.weir.fst", 
                  sep = "\t", header=TRUE)

##########################################################
####Generate a function to impute NA in our .fst files####
##########################################################

process_dataframe <- function(data) {
  # Get the unique chromosomes
  chromosomes <- unique(data$CHROM)
  
  # Initialize an empty list to store processed data for each chromosome
  processed_data_list <- list()
  
  # Loop through each chromosome
  for (chrom in chromosomes) {
    # Subset the data for the current chromosome
    chrom_data <- data %>% filter(CHROM == chrom)
    
    # Create a complete sequence of BIN_START values for the current chromosome
    min_bin <- min(chrom_data$BIN_START)
    max_bin <- max(chrom_data$BIN_START)
    bin_step <- 2500 # Assuming BIN_START increments by 2500
    
    min_bin_end <- min(chrom_data$BIN_END)
    max_bin_end <- max(chrom_data$BIN_END)
    bin_step_end <- 2500 # Assuming BIN_START increments by 2500
    
    # Generate a complete sequence of BIN_START values
    complete_bins <- seq(min_bin, max_bin, by = bin_step)
    complete_bins_end <- seq(min_bin_end, max_bin_end, by = bin_step_end)
    
    # Create a data frame with the complete sequence
    complete_data <- data.frame(BIN_START = complete_bins, BIN_END = complete_bins_end)
    
    # Merge the complete sequence with the chromosome data
    merged_data <- complete_data %>%
      left_join(chrom_data, by = c("BIN_START", "BIN_END"))
    
    # Fill the missing CHROM values with the current chromosome
    merged_data$CHROM <- chrom
    
    # Impute missing FST values with 0
    merged_data <- merged_data %>%
      mutate(WEIGHTED_FST = ifelse(is.na(WEIGHTED_FST), 0, WEIGHTED_FST),
             MEAN_FST = ifelse(is.na(MEAN_FST), 0, MEAN_FST),
             N_VARIANTS = ifelse(is.na(N_VARIANTS), 0, N_VARIANTS))
    
    # Append the processed data for the current chromosome to the list
    processed_data_list[[chrom]] <- merged_data
  }
  
  # Combine the processed data for all chromosomes
  processed_data <- bind_rows(processed_data_list)
  
  return(processed_data)
}


file_names <- c("df1pop", "df1chr2", "df1chr4", "df1chr11P", "df1chr16")

# Loop through each data frame
for (df_name in file_names) {
  # Get the data frame from its name
  data <- get(df_name)
  
  # Process and impute the missing data
  imputed_data <- process_dataframe(data)
  
  assign(df_name, imputed_data)
  
}

rm(imputed_data, data)

df1pop <- df1pop[,c(3,1,5)]
df1pop <- cbind(df1pop, "SNP"= 1:nrow(df1pop)) 
df1pop$CHROM<-gsub("chrPGA_chromosome_","",as.character(df1pop$CHROM))
df1pop <- df1pop %>% select(SNP, everything())
df1pop$SNP <- as.character(df1pop$SNP)
colnames(df1pop) <- c("SNP", "chr", "position", "globalFST")

df1chr2 <- df1chr2[,c(3,1,5)]
df1chr2 <- cbind(df1chr2, "SNP"= 1:nrow(df1chr2)) 
df1chr2$CHROM<-gsub("chrPGA_chromosome_","",as.character(df1chr2$CHROM))
df1chr2 <- df1chr2 %>% select(SNP, everything())
df1chr2$SNP <- as.character(df1chr2$SNP)
colnames(df1chr2) <- c("SNP", "chr", "position", "globalFST")

df1chr4 <- df1chr4[,c(3,1,5)]
df1chr4 <- cbind(df1chr4, "SNP"= 1:nrow(df1chr4)) 
df1chr4$CHROM<-gsub("chrPGA_chromosome_","",as.character(df1chr4$CHROM))
df1chr4 <- df1chr4 %>% select(SNP, everything())
df1chr4$SNP <- as.character(df1chr4$SNP)
colnames(df1chr4) <- c("SNP", "chr", "position", "globalFST")

df1chr11P <- df1chr11P[,c(3,1,5)]
df1chr11P <- cbind(df1chr11P, "SNP"= 1:nrow(df1chr11P)) 
df1chr11P$CHROM<-gsub("chrPGA_chromosome_","",as.character(df1chr11P$CHROM))
df1chr11P <- df1chr11P %>% select(SNP, everything())
df1chr11P$SNP <- as.character(df1chr11P$SNP)
colnames(df1chr11P) <- c("SNP", "chr", "position", "globalFST")

df1chr16 <- df1chr16[,c(3,1,5)]
df1chr16 <- cbind(df1chr16, "SNP"= 1:nrow(df1chr16)) 
df1chr16$CHROM<-gsub("chrPGA_chromosome_","",as.character(df1chr16$CHROM))
df1chr16 <- df1chr16 %>% select(SNP, everything())
df1chr16$SNP <- as.character(df1chr16$SNP)
colnames(df1chr16) <- c("SNP", "chr", "position", "globalFST")


#### Find all diff and corr genes cooridnates####
overlap <- intersect(Corr$`Corr$Row`, Diff$V1)

Diff <- Diff %>% filter(!(V1 %in% overlap))
Corr_gff <- inner_join(gff, Corr, by = c("gene" = "Corr$Row"))
Diff_gff <- inner_join(gff, Diff, by = c("gene" = "V1"))

Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]

write.table(Corr_gff, "gff_corr_TOT.txt")
write.table(Diff_gff, "gff_diff_TOT.txt")


DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1pop[df1pop$chr == Diff_gff$Chrom[i] & df1pop$position >= Diff_gff$Start[i]-1250 & df1pop$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- as.vector(c(t(DiffGene)))


CorrGene <- data.frame()
for (i in 1:nrow(Corr_gff)){
  sweep <- as.data.frame(df1pop[df1pop$chr == Corr_gff$Chrom[i] & df1pop$position >= Corr_gff$Start[i]-1250 & df1pop$position <= Corr_gff$Finish[i]+1250,1])
  CorrGene <- rbind.fill(CorrGene,sweep)
}
CorrGene <- as.vector(c(t(CorrGene)))


#### FST manhattan plot all ####
CMplot(df1pop, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene, CorrGene), 
       highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
       highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("FST"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=4, width=10, file.name="GlobalFST_allgenes_Pop")



#### FST violin plots all ####
df1pop$distr <- c("Genomic")
Diff_dens <- df1pop %>% filter((SNP %in% DiffGene))
Corr_dens <- df1pop %>% filter((SNP %in% CorrGene))

pdf(file="FST_globalpop_allwindows.pdf", height=4, width=3)
ggplot(df1pop, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1)+
  scale_color_manual(values="grey80", "black")+
  ylim(c(-0.2,1))+
  theme_classic()
dev.off()

pdf(file="FST_globalpop_DiffWindows.pdf", height=4, width=3)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1)+
  scale_color_manual(values=c("#FFBF00", "black"))+
  ylim(c(-0.2,1))+
  theme_classic()
dev.off()

pdf(file="FST_globalpop_CorrWindows.pdf", height=4, width=3)
ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#CC5500")+
  geom_boxplot(color = "black", width = 0.1)+
  scale_color_manual(values="#CC5500", "black")+
  theme(legend.position = "none")+
  ylim(c(-0.2,1))+
  theme_classic()
dev.off()

all_allvsDiff_ttest <- wilcox.test( x = df1pop$globalFST, y = Diff_dens$globalFST, alternative = "two.sided", paired=F)
all_allvsCorr_ttest <- wilcox.test( x = df1pop$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
all_DiffvsCorr_ttest <- wilcox.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")

all_allvsDiff_kstest <- ks.test( x = df1pop$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
all_allvsCorr_kstest <- ks.test( x = df1pop$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
all_DiffvsCorr_kstest <- ks.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")



all_allvsDiff_ttest <- data.frame(statistic = all_allvsDiff_ttest$statistic, p_value = all_allvsDiff_ttest$p.value, df = all_allvsDiff_ttest$parameter, method = all_allvsDiff_ttest$method)
all_allvsCorr_ttest <- data.frame(statistic = all_allvsCorr_ttest$statistic, p_value = all_allvsCorr_ttest$p.value, df = all_allvsCorr_ttest$parameter, method = all_allvsCorr_ttest$method)
all_DiffvsCorr_ttest <- data.frame(statistic = all_DiffvsCorr_ttest$statistic, p_value = all_DiffvsCorr_ttest$p.value, df = all_DiffvsCorr_ttest$parameter, method = all_DiffvsCorr_ttest$method)

all_allvsDiff_kstest <- data.frame(statistic = all_allvsDiff_kstest$statistic, p_value = all_allvsDiff_kstest$p.value,  method = all_allvsDiff_kstest$method)
all_allvsCorr_kstest <- data.frame(statistic = all_allvsCorr_kstest$statistic, p_value = all_allvsCorr_kstest$p.value, method = all_allvsCorr_kstest$method)
all_DiffvsCorr_kstest <- data.frame(statistic = all_DiffvsCorr_kstest$statistic, p_value = all_DiffvsCorr_kstest$p.value, method = all_DiffvsCorr_kstest$method)


write.table(all_allvsDiff_ttest, "all_allvsDiff_ttest.txt", row.names = FALSE)
write.table(all_allvsCorr_ttest, "all_allvsCorr_ttest.txt", row.names = FALSE)
write.table(all_DiffvsCorr_ttest, "all_DiffvsCorr_ttest.txt", row.names = FALSE)

write.table(all_allvsDiff_kstest, "all_allvsDiff_kstest.txt", row.names = FALSE)
write.table(all_allvsCorr_kstest, "all_allvsCorr_kstest.txt", row.names = FALSE)
write.table(all_DiffvsCorr_kstest, "all_DiffvsCorr_kstest.txt", row.names = FALSE)


##### Find Pop diff and corr genes cooridnates ####
Cpop <- as.data.frame(Cpop$Row)
overlap <- intersect(Cpop$`Cpop$Row`, Dpop$V1)
Dpop <- Dpop %>% filter(!(V1 %in% overlap))
Corr_gff <- inner_join(gff, Cpop, by = c("gene" = 'Cpop$Row'))
Diff_gff <- inner_join(gff, Dpop, by = c("gene" = "V1"))

Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]

write.table(Corr_gff, "gff_corr_pop.txt")
write.table(Diff_gff, "gff_diff_pop.txt")

DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1pop[df1pop$chr == Diff_gff$Chrom[i] & df1pop$position >= Diff_gff$Start[i]-1250 & df1pop$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- c(t(DiffGene))

CorrGene <- data.frame()
for (n in 1:nrow(Corr_gff)){
  sweep <- as.data.frame(df1pop[df1pop$chr == Corr_gff$Chrom[n] & df1pop$position >= Corr_gff$Start[n]-1250 & df1pop$position <= Corr_gff$Finish[n]+1250,1])
  CorrGene <- rbind.fill(CorrGene,sweep)
}
CorrGene <- c(t(CorrGene))

df1pop <- df1pop[,1:4]

##### FST manhattan plot pop ####
CMplot(df1pop, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene, CorrGene), 
       highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
       highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("Pop_genes"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=3, width=10, file.name="GlobalFST_Popgenes")

##### FST violin plots pop ####
df1pop$distr <- c("Genomic")
Diff_dens <- df1pop %>% filter((SNP %in% DiffGene))
Corr_dens <- df1pop %>% filter((SNP %in% CorrGene))

pdf(file="FST_pop_allwindows.pdf", height=2, width=2.4)
ggplot(df1pop, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="grey80", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_pop_DiffWindows.pdf", height=2, width=2.4)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#FFBF00", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_pop_CorrWindows.pdf", height=2, width=2.4)
ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#CC5500")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#CC5500", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pop_allvsDiff_ttest <- wilcox.test( x = df1pop$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
pop_allvsCorr_ttest <- wilcox.test( x = df1pop$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
pop_DiffvsCorr_ttest <- wilcox.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")

pop_allvsDiff_kstest <- ks.test( x = df1pop$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
pop_allvsCorr_kstest <- ks.test( x = df1pop$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
pop_DiffvsCorr_kstest <- ks.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")


pop_allvsDiff_ttest <- data.frame(statistic = pop_allvsDiff_ttest$statistic, p_value = pop_allvsDiff_ttest$p.value, df = pop_allvsDiff_ttest$parameter, method = pop_allvsDiff_ttest$method)
pop_allvsCorr_ttest <- data.frame(statistic = pop_allvsCorr_ttest$statistic, p_value = pop_allvsCorr_ttest$p.value, df = pop_allvsCorr_ttest$parameter, method = pop_allvsCorr_ttest$method)
pop_DiffvsCorr_ttest <- data.frame(statistic = pop_DiffvsCorr_ttest$statistic, p_value = pop_DiffvsCorr_ttest$p.value, df = pop_DiffvsCorr_ttest$parameter, method = pop_DiffvsCorr_ttest$method)

pop_allvsDiff_kstest <- data.frame(statistic = pop_allvsDiff_kstest$statistic, p_value = pop_allvsDiff_kstest$p.value, method = pop_allvsDiff_kstest$method)
pop_allvsCorr_kstest <- data.frame(statistic = pop_allvsCorr_kstest$statistic, p_value = pop_allvsCorr_kstest$p.value, method = pop_allvsCorr_kstest$method)
pop_DiffvsCorr_kstest <- data.frame(statistic = pop_DiffvsCorr_kstest$statistic, p_value = pop_DiffvsCorr_kstest$p.value, method = pop_DiffvsCorr_kstest$method)


write.table(pop_allvsDiff_ttest, "pop_allvsDiff_ttest.txt", row.names = FALSE)
write.table(pop_allvsCorr_ttest, "pop_allvsCorr_ttest.txt", row.names = FALSE)
write.table(pop_DiffvsCorr_ttest, "pop_DiffvsCorr_ttest.txt", row.names = FALSE)

write.table(pop_allvsDiff_kstest, "pop_allvsDiff_kstest.txt", row.names = FALSE)
write.table(pop_allvsCorr_kstest, "pop_allvsCorr_kstest.txt", row.names = FALSE)
write.table(pop_DiffvsCorr_kstest, "pop_DiffvsCorr_kstest.txt", row.names = FALSE)



##### Find Chr2 diff and corr genes cooridnates ####
Cchr2 <- as.data.frame(Cchr2$Row)
overlap <- intersect(Cchr2$`Cchr2$Row`, Dchr2$V1)
Dchr2 <- Dchr2 %>% filter(!(V1 %in% overlap))
Corr_gff <- inner_join(gff, Cchr2, by = c("gene" = "Cchr2$Row"))
Diff_gff <- inner_join(gff, Dchr2, by = c("gene" = "V1"))

Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]

DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1chr2[df1chr2$chr == Diff_gff$Chrom[i] & df1chr2$position >= Diff_gff$Start[i]-1250 & df1chr2$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- c(t(DiffGene))

CorrGene <- data.frame()
for (n in 1:nrow(Corr_gff)){
  nsweep <- as.data.frame(df1chr2[df1chr2$chr == Corr_gff$Chrom[n] & df1chr2$position >= Corr_gff$Start[n]-1250 & df1chr2$position <= Corr_gff$Finish[n]+1250,1])
  CorrGene <- rbind.fill(CorrGene,nsweep)
}
CorrGene <- c(t(CorrGene))

write.table(Corr_gff, "gff_corr_chr2.txt")
write.table(Diff_gff, "gff_diff_chr2.txt")


df1chr2 <- df1chr2[,1:4]
##### FST manhattan plot chr2 ####
CMplot(df1chr2, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene, CorrGene), 
       highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
       highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("Pop_genes"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=3, width=10, file.name="Chr2FST_chr2genes")

##### FST violin plots chr2 ####
df1chr2$distr <- c("Genomic")
Diff_dens <- df1chr2 %>% filter((SNP %in% DiffGene))
Corr_dens <- df1chr2 %>% filter((SNP %in% CorrGene))

pdf(file="FST_chr2_allwindows.pdf", height=2, width=2.4)
ggplot(df1chr2, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="grey80", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr2_DiffWindows.pdf", height=2, width=2.4)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#FFBF00", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr2_CorrWindows.pdf", height=2, width=2.4)
ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#CC5500")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#CC5500", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

chr2_allvsDiff_ttest <- wilcox.test( x = df1chr2$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr2_allvsCorr_ttest <- wilcox.test( x = df1chr2$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr2_DiffvsCorr_ttest <- wilcox.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")

chr2_allvsDiff_kstest <- ks.test( x = df1chr2$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr2_allvsCorr_kstest <- ks.test( x = df1chr2$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr2_DiffvsCorr_kstest <- ks.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")


chr2_allvsDiff_ttest <- data.frame(statistic = chr2_allvsDiff_ttest$statistic, p_value = chr2_allvsDiff_ttest$p.value, df = chr2_allvsDiff_ttest$parameter, method = chr2_allvsDiff_ttest$method)
chr2_allvsCorr_ttest <- data.frame(statistic = chr2_allvsCorr_ttest$statistic, p_value = chr2_allvsCorr_ttest$p.value, df = chr2_allvsCorr_ttest$parameter, method = chr2_allvsCorr_ttest$method)
chr2_DiffvsCorr_ttest <- data.frame(statistic = chr2_DiffvsCorr_ttest$statistic, p_value = chr2_DiffvsCorr_ttest$p.value, df = chr2_DiffvsCorr_ttest$parameter, method = chr2_DiffvsCorr_ttest$method)

chr2_allvsDiff_kstest <- data.frame(statistic = chr2_allvsDiff_kstest$statistic, p_value = chr2_allvsDiff_kstest$p.value, method = chr2_allvsDiff_kstest$method)
chr2_allvsCorr_kstest <- data.frame(statistic = chr2_allvsCorr_kstest$statistic, p_value = chr2_allvsCorr_kstest$p.value, method = chr2_allvsCorr_kstest$method)
chr2_DiffvsCorr_kstest <- data.frame(statistic = chr2_DiffvsCorr_kstest$statistic, p_value = chr2_DiffvsCorr_kstest$p.value, method = chr2_DiffvsCorr_kstest$method)


write.table(chr2_allvsDiff_ttest, "chr2_allvsDiff_ttest.txt", row.names = FALSE)
write.table(chr2_allvsCorr_ttest, "chr2_allvsCorr_ttest.txt", row.names = FALSE)
write.table(chr2_DiffvsCorr_ttest, "chr2_DiffvsCorr_ttest.txt", row.names = FALSE)

write.table(chr2_allvsDiff_kstest, "chr2_allvsDiff_kstest.txt", row.names = FALSE)
write.table(chr2_allvsCorr_kstest, "chr2_allvsCorr_kstest.txt", row.names = FALSE)
write.table(chr2_DiffvsCorr_kstest, "chr2_DiffvsCorr_kstest.txt", row.names = FALSE)


##### Find Chr4 diff and corr genes cooridnates ####
Cchr4 <- as.data.frame(Cchr4$Row)
overlap <- intersect(Cchr4$`Cchr4$Row`, Dchr4$V1)
Dchr4 <- Dchr4 %>% filter(!(V1 %in% overlap))
Corr_gff <- inner_join(gff, Cchr4, by = c("gene" = "Cchr4$Row"))
Diff_gff <- inner_join(gff, Dchr4, by = c("gene" = "V1"))

Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]


write.table(Corr_gff, "gff_corr_chr4.txt")
write.table(Diff_gff, "gff_diff_chr4.txt")


DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1chr4[df1chr4$chr == Diff_gff$Chrom[i] & df1chr4$position >= Diff_gff$Start[i]-1250 & df1chr4$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- c(t(DiffGene))

CorrGene <- data.frame()
for (n in 1:nrow(Corr_gff)){
  nsweep <- as.data.frame(df1chr4[df1chr4$chr == Corr_gff$Chrom[n] & df1chr4$position >= Corr_gff$Start[n]-1250 & df1chr4$position <= Corr_gff$Finish[n]+1250,1])
  CorrGene <- rbind.fill(CorrGene,nsweep)
}
CorrGene <- c(t(CorrGene))

##### FST manhattan plot chr4 ####
CMplot(df1chr4, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene, CorrGene), 
       highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
       highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("Pop_genes"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=3, width=10, file.name="Chr4FST_chr4genes")

##### FST violin plots chr4 ####
df1chr4$distr <- c("Genomic")
Diff_dens <- df1chr4 %>% filter((SNP %in% DiffGene))
Corr_dens <- df1chr4 %>% filter((SNP %in% CorrGene))

pdf(file="FST_chr4_allwindows.pdf", height=2, width=2.4)
ggplot(df1chr4, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="grey80", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr4_DiffWindows.pdf", height=2, width=2.4)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#FFBF00", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr4_CorrWindows.pdf", height=2, width=2.4)
ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#CC5500")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#CC5500", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

chr4_allvsDiff_ttest <- wilcox.test( x = df1chr4$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr4_allvsCorr_ttest <- wilcox.test( x = df1chr4$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr4_DiffvsCorr_ttest <- wilcox.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")

chr4_allvsDiff_kstest <- ks.test( x = df1chr4$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr4_allvsCorr_kstest <- ks.test( x = df1chr4$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr4_DiffvsCorr_kstest <- ks.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")



chr4_allvsDiff_ttest <- data.frame(statistic = chr4_allvsDiff_ttest$statistic, p_value = chr4_allvsDiff_ttest$p.value, df = chr4_allvsDiff_ttest$parameter, method = chr4_allvsDiff_ttest$method)
chr4_allvsCorr_ttest <- data.frame(statistic = chr4_allvsCorr_ttest$statistic, p_value = chr4_allvsCorr_ttest$p.value, df = chr4_allvsCorr_ttest$parameter, method = chr4_allvsCorr_ttest$method)
chr4_DiffvsCorr_ttest <- data.frame(statistic = chr4_DiffvsCorr_ttest$statistic, p_value = chr4_DiffvsCorr_ttest$p.value, df = chr4_DiffvsCorr_ttest$parameter, method = chr4_DiffvsCorr_ttest$method)

chr4_allvsDiff_kstest <- data.frame(statistic = chr4_allvsDiff_kstest$statistic, p_value = chr4_allvsDiff_kstest$p.value, method = chr4_allvsDiff_kstest$method)
chr4_allvsCorr_kstest <- data.frame(statistic = chr4_allvsCorr_kstest$statistic, p_value = chr4_allvsCorr_kstest$p.value,  method = chr4_allvsCorr_kstest$method)
chr4_DiffvsCorr_kstest <- data.frame(statistic = chr4_DiffvsCorr_kstest$statistic, p_value = chr4_DiffvsCorr_kstest$p.value, method = chr4_DiffvsCorr_kstest$method)


write.table(chr4_allvsDiff_ttest, "chr4_allvsDiff_ttest.txt", row.names = FALSE)
write.table(chr4_allvsCorr_ttest, "chr4_allvsCorr_ttest.txt", row.names = FALSE)
write.table(chr4_DiffvsCorr_ttest, "chr4_DiffvsCorr_ttest.txt", row.names = FALSE)

write.table(chr4_allvsDiff_kstest, "chr4_allvsDiff_kstest.txt", row.names = FALSE)
write.table(chr4_allvsCorr_kstest, "chr4_allvsCorr_kstest.txt", row.names = FALSE)
write.table(chr4_DiffvsCorr_kstest, "chr4_DiffvsCorr_kstest.txt", row.names = FALSE)



##### Find Chr11 diff and corr genes cooridnates ####
Cchr11 <- as.data.frame(Cchr11$Row)
overlap <- intersect(Cchr11$`Cchr11$Row`, Dchr11$V1)
Dchr11 <- Dchr11 %>% filter(!(V1 %in% overlap))
Corr_gff <- inner_join(gff, Cchr11, by = c("gene" = "Cchr11$Row"))
Diff_gff <- inner_join(gff, Dchr11, by = c("gene" = "V1"))

Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]

write.table(Corr_gff, "gff_corr_chr11.txt")
write.table(Diff_gff, "gff_diff_chr11.txt")

# DiffGene <- data.frame()
# for (i in 1:nrow(Diff_gff)){
#   sweep <- as.data.frame(df1chr11G[df1chr11G$chr == Diff_gff$Chrom[i] & df1chr11G$position >= Diff_gff$Start[i]-1250 & df1chr11G$position <= Diff_gff$Finish[i]+1250,1])
#   DiffGene <- rbind.fill(DiffGene,sweep)
# }
# DiffGene <- c(t(DiffGene))
# 
# CorrGene <- data.frame()
# for (n in 1:nrow(Corr_gff)){
#   nsweep <- as.data.frame(df1chr11G[df1chr11G$chr == Corr_gff$Chrom[n] & df1chr11G$position >= Corr_gff$Start[n]-1250 & df1chr11G$position <= Corr_gff$Finish[n]+1250,1])
#   CorrGene <- rbind.fill(CorrGene,nsweep)
# }
# CorrGene <- c(t(CorrGene))

##### FST manhattan plot chr11 ####
# CMplot(df1chr11G, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
#        #threshold=list(c(0.73)),
#        ylim=c(-0.25,1),
#        signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
#        highlight=c(DiffGene, CorrGene), 
#        highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
#        highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
#        # threshold.col=c("yellow1"), signal.col=c("#d35400"),
#        chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
#        ylab=c("Pop_genes"),
#        chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
#        file="pdf", height=4, width=10, file.name="Chr11FST_chr11Ggenes")
# 
# ##### FST violin plots chr11 ####
# df1chr11G$distr <- c("Genomic")
# Diff_dens <- df1chr11G %>% filter((SNP %in% DiffGene))
# Corr_dens <- df1chr11G %>% filter((SNP %in% CorrGene))
# 
# pdf(file="FST_chr11G_allwindows.pdf", height=2, width=2.4)
# ggplot(df1chr11G, aes(x=distr, y= globalFST, fill=distr, color=distr))+
#   geom_violin(alpha=0.5, lwd=0.5)+
#   scale_fill_manual(values="grey80")+
#   geom_boxplot(color = "black", width = 0.1)+
#   scale_color_manual(values="grey80", "black")+
#   ylim(c(-0.2,1))+
#   theme_classic()
# dev.off()
# 
# pdf(file="FST_chr11G_DiffWindows.pdf", height=2, width=2.4)
# ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
#   geom_violin(alpha=0.5, lwd=0.5)+
#   scale_fill_manual(values="#FFBF00")+
#   geom_boxplot(color = "black", width = 0.1)+
#   scale_color_manual(values=c("#FFBF00", "black"))+
#   ylim(c(-0.2,1))+
#   theme_classic()
# dev.off()
# 
# pdf(file="FST_chr11G_CorrWindows.pdf", height=2, width=2.4)
# ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
#   geom_violin(alpha=0.5, lwd=0.5)+
#   scale_fill_manual(values="#CC5500")+
#   geom_boxplot(color = "black", width = 0.1)+
#   scale_color_manual(values="#CC5500", "black")+
#   ylim(c(-0.2,1))+
#   theme_classic()
# dev.off()



##### Find Chr11 diff and corr genes cooridnates ####
DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1chr11P[df1chr11P$chr == Diff_gff$Chrom[i] & df1chr11P$position >= Diff_gff$Start[i]-1250 & df1chr11P$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- c(t(DiffGene))

CorrGene <- data.frame()
for (n in 1:nrow(Corr_gff)){
  nsweep <- as.data.frame(df1chr11P[df1chr11P$chr == Corr_gff$Chrom[n] & df1chr11P$position >= Corr_gff$Start[n]-1250 & df1chr11P$position <= Corr_gff$Finish[n]+1250,1])
  CorrGene <- rbind.fill(CorrGene,nsweep)
}
CorrGene <- c(t(CorrGene))

##### FST manhattan plot chr11 ####
CMplot(df1chr11P, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene, CorrGene), 
       highlight.col=c(rep("#FFBF00", times = length(DiffGene)), rep("#CC5500",times=length(CorrGene))),
       highlight.cex = c(0.7), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("Pop_genes"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=3, width=10, file.name="Chr11FST_chr11Pgenes")

##### FST violin plots chr11 ####
df1chr11P$distr <- c("Genomic")
Diff_dens <- df1chr11P %>% filter((SNP %in% DiffGene))
Corr_dens <- df1chr11P %>% filter((SNP %in% CorrGene))

pdf(file="FST_chr11P_allwindows.pdf", height=2, width=2.4)
ggplot(df1chr11P, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="grey80", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr11P_DiffWindows.pdf", height=2, width=2.4)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#FFBF00", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr11P_CorrWindows.pdf", height=2, width=2.4)
ggplot(Corr_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#CC5500")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#CC5500", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

shapiro.test(df1chr11P$globalFST)   
shapiro.test(Diff_dens$globalFST)   
shapiro.test(Corr_dens$globalFST)   

chr11_allvsDiff_ttest <- wilcox.test( x = df1chr11P$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr11_allvsCorr_ttest <- wilcox.test( x = df1chr11P$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr11_DiffvsCorr_ttest <- wilcox.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")

chr11_allvsDiff_kstest <- ks.test( x = df1chr11P$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")
chr11_allvsCorr_kstest <- ks.test( x = df1chr11P$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")
chr11_DiffvsCorr_kstest <- ks.test( x = Diff_dens$globalFST, y = Corr_dens$globalFST, alternative = "two.sided")



chr11_allvsDiff_ttest <- data.frame(statistic = chr11_allvsDiff_ttest$statistic, p_value = chr11_allvsDiff_ttest$p.value, df = chr11_allvsDiff_ttest$parameter, method = chr11_allvsDiff_ttest$method)
chr11_allvsCorr_ttest <- data.frame(statistic = chr11_allvsCorr_ttest$statistic, p_value = chr11_allvsCorr_ttest$p.value, df = chr11_allvsCorr_ttest$parameter, method = chr11_allvsCorr_ttest$method)
chr11_DiffvsCorr_ttest <- data.frame(statistic = chr11_DiffvsCorr_ttest$statistic, p_value = chr11_DiffvsCorr_ttest$p.value, df = chr11_DiffvsCorr_ttest$parameter, method = chr11_DiffvsCorr_ttest$method)

chr11_allvsDiff_kstest <- data.frame(statistic = chr11_allvsDiff_kstest$statistic, p_value = chr11_allvsDiff_kstest$p.value, method = chr11_allvsDiff_kstest$method)
chr11_allvsCorr_kstest <- data.frame(statistic = chr11_allvsCorr_kstest$statistic, p_value = chr11_allvsCorr_kstest$p.value, method = chr11_allvsCorr_kstest$method)
chr11_DiffvsCorr_kstest <- data.frame(statistic = chr11_DiffvsCorr_kstest$statistic, p_value = chr11_DiffvsCorr_kstest$p.value, method = chr11_DiffvsCorr_kstest$method)


write.table(chr11_allvsDiff_ttest, "chr11_allvsDiff_ttest.txt", row.names = FALSE)
write.table(chr11_allvsCorr_ttest, "chr11_allvsCorr_ttest.txt", row.names = FALSE)
write.table(chr11_DiffvsCorr_ttest, "chr11_DiffvsCorr_ttest.txt", row.names = FALSE)

write.table(chr11_allvsDiff_kstest, "chr11_allvsDiff_kstest.txt", row.names = FALSE)
write.table(chr11_allvsCorr_kstest, "chr11_allvsCorr_kstest.txt", row.names = FALSE)
write.table(chr11_DiffvsCorr_kstest, "chr11_DiffvsCorr_kstest.txt", row.names = FALSE)

##### Find Chr16 diff and corr genes cooridnates ####
# Cchr16 <- as.data.frame(Cchr16$Row)
# overlap <- intersect(Cchr16$`Cchr16$Row`, Dchr16$V1)
# Dchr16 <- Dchr16 %>% filter(!(V1 %in% overlap))
# Corr_gff <- inner_join(gff, Cchr16, by = c("gene" = "Cchr16$Row"))
Diff_gff <- inner_join(gff, Dchr16, by = c("gene" = "V1"))

#Corr_gff <- Corr_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))
Diff_gff <- Diff_gff %>% group_by(gene) %>%  dplyr::summarise( Start = min(Start), Finish = max(Finish))

#Corr_chr<- unique(gff[gff$gene %in% Corr_gff$gene, c(1,4)])
Diff_chr<- unique(gff[gff$gene %in% Diff_gff$gene, c(1,4)])

#Corr_gff <- as.data.frame(inner_join(Corr_gff, Corr_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]
Diff_gff <- as.data.frame(inner_join(Diff_gff, Diff_chr, by = c("gene" = "gene")))[,c(4,2,3,1)]

write.table(Diff_gff, "gff_diff_chr16.txt")

DiffGene <- data.frame()
for (i in 1:nrow(Diff_gff)){
  sweep <- as.data.frame(df1chr16[df1chr16$chr == Diff_gff$Chrom[i] & df1chr16$position >= Diff_gff$Start[i]-1250 & df1chr16$position <= Diff_gff$Finish[i]+1250,1])
  DiffGene <- rbind.fill(DiffGene,sweep)
}
DiffGene <- c(t(DiffGene))

# CorrGene <- data.frame()
# for (n in 1:nrow(Corr_gff)){
#   nsweep <- as.data.frame(df1chr16[df1chr16$chr == Corr_gff$Chrom[n] & df1chr16$position >= Corr_gff$Start[n] & df1chr16$position <= Corr_gff$Finish[n],1]
#   CorrGene <- rbind.fill(DiffGene,nsweep)
# }
# CorrGene <- c(t(CorrGene))

##### FST manhattan plot chr16 ####
CMplot(df1chr16, plot.type="m", type="p",r=0.5, band=1, outward=T, col=c("grey80","grey80"),
       #threshold=list(c(0.73)),
       ylim=c(-0.25,1),
       signal.line = NULL,  cir.band = 0.1, LOG10=F, chr.border=FALSE,
       highlight=c(DiffGene), 
       highlight.col=c("#FFBF00"),
       highlight.cex = c(0.5), highlight.type = "p", #highlight.pch = c(15),
       # threshold.col=c("yellow1"), signal.col=c("#d35400"),
       chr.labels=paste("Chr",c(1:16),sep=""), #signal.cex=c(0.35),
       ylab=c("Pop_genes"),
       chr.den.col=c("lightblue","lightblue", "#1f618d", "black"), cex=c(0.2), points.alpha=(100),
       file="pdf", height=3, width=10, file.name="Chr16FST_chr16genes")

##### FST violin plots chr16 ####
df1chr16$distr <- c("Genomic")
Diff_dens <- df1chr16 %>% filter((SNP %in% DiffGene))

pdf(file="FST_chr16_allwindows.pdf", height=2, width=2.4)
ggplot(df1chr16, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="grey80")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="grey80", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()

pdf(file="FST_chr16_DiffWindows.pdf", height=2, width=2.4)
ggplot(Diff_dens, aes(x=distr, y= globalFST, fill=distr, color=distr))+
  geom_violin(alpha=0.5, lwd=0.5)+
  scale_fill_manual(values="#FFBF00")+
  geom_boxplot(color = "black", width = 0.1, outlier.shape = NA )+
  scale_color_manual(values="#FFBF00", "black")+
  scale_y_continuous(limits = c(-0.25, 1), breaks = seq(-0.2, 1, by = 0.2)) +
  theme_classic()+
  theme(legend.position="none")
dev.off()


chr16_allvsDiff_ttest <- wilcox.test( x = df1chr16$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")

chr16_allvsDiff_kstest <- ks.test( x = df1chr16$globalFST, y = Diff_dens$globalFST, alternative = "two.sided")


chr16_allvsDiff_ttest <- data.frame(statistic = chr16_allvsDiff_ttest$statistic, p_value = chr16_allvsDiff_ttest$p.value, df = chr16_allvsDiff_ttest$parameter, method = chr16_allvsDiff_ttest$method)

chr16_allvsDiff_kstest <- data.frame(statistic = chr16_allvsDiff_kstest$statistic, p_value = chr16_allvsDiff_kstest$p.value,  method = chr16_allvsDiff_kstest$method)

write.table(chr16_allvsDiff_ttest, "chr16_allvsDiff_ttest.txt", row.names = FALSE)

write.table(chr16_allvsDiff_kstest, "chr16_allvsDiff_kstest.txt", row.names = FALSE)

