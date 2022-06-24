

library(dplyr)
library(readr)
library(tidyr)
library(sqldf)
library(stringr)

##* path
# RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/"
# RES="Z:/projects/Matei_lab/chip_cut_run/data/res/"
RES="/Users/pallavi/Documents/k/projects/Matei_lab/chip_cut_run/data/res/diffbind_merged/annotated_peaks/"
path = "/Users/pallavi/Documents/k/projects/Matei_lab/chip_cut_run/data/res/diffbind_merged/alternate_promoters/"


all = read_tsv(paste0(RES, "k4me3.csv"))
all = subset(all, annotation=="Promoter" | annotation=="Distal Intergenic")

all_sub = all[, c(1:3,19,21:22, 9:11)]
all_sub = rename(all_sub, chr = seqnames)

# --------------------------------------------------------
#find alternate promoters
# --------------------------------------------------------

#remove na's of genes
all_sub = subset(all_sub, !is.na(ENSEMBL))

# #pvalue cutoff - already considered before
# all_sub = subset(all_sub, p.value < 0.01)

# fold change is +
pos = subset(all_sub, Fold > 0)
pos1 = pos[order(pos$ENSEMBL, decreasing = FALSE),]
pos1$gene_tc <- paste(pos1$ENSEMBL, pos1$transcriptId, sep="_")


# fold change is -
neg = subset(all_sub, Fold < 0)
neg1 = neg[order(neg$ENSEMBL, decreasing = FALSE),]
neg1$gene_tc <- paste(neg1$ENSEMBL, neg1$transcriptId, sep="_")

#- alternate promoters
me <- merge(pos1, neg1, by.x = c("ENSEMBL"), by.y=c("ENSEMBL"), all = FALSE)
me$SYMBOL.y = NULL
write.table(me, file = paste0(path, "alt_promoter", "k4me3.csv"), sep="\t", quote=F, row.names=F)

