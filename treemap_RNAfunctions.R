#install.packages("treemap")
library(treemap)

setwd("E:/Splicata/16S-RNA-DNA/RNA_function")

treemap_MF_Diff <- read.table("Revigo_MF_Corr_onlychr.tsv", sep="\t", header=T)
treemap_BP_Diff <- read.table("Revigo_BP_Corr_onlychr.tsv", sep="\t", header=T)

header<- c("GO term",	"Description",	"Functional Category",	"LogSize",	"Frequency",	"Uniqueness", "Dispensability",	"Representative")

colnames(treemap_MF_Diff) <- header
colnames(treemap_BP_Diff) <- header

pdf("treemap_MF_Corr_onlychr.pdf",width=10.5, height=4.3)
treemap(treemap_MF_Diff, index=c("Functional Category","Description"), vSize="Uniqueness", 
        fontsize.labels=c(10,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("black","grey50"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=T,                        # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","grey30"),             # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(3,1),                        # Width of colors
        palette = "PuRd",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="",                      # Customize your title
        #fontsize.title=12,                       # Size of the title
        )
dev.off()


pdf("treemap_BP_Corr_onlychr.pdf",width=8.5, height=4.3)
treemap(treemap_BP_Diff, index=c("Functional Category","Description"), vSize="Uniqueness", type = "index",
        fontsize.labels=c(10,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("black","grey50"),    # Color of labels
        fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("center", "center"), 
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=T,                        # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","grey30"),             # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(3,1),                        # Width of colors
        palette = "-Greens",                        # Select your color palette from the RColorBrewer presets or make your own.
        title="",                      # Customize your title
        #fontsize.title=12,                       # Size of the title
)
dev.off()
