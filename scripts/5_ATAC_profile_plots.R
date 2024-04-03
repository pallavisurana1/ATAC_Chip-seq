
#- conda activate R411
##* packages
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(ggplot2)
library(GenomicRanges)
# library(ReactomePA)
library(BRGenomics)


###* path
RES="/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/"

setwd(RES)

# ##* get files
# file_bed = "/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/overlap_peaks/bed_FTO_PCDH_overlap.bed"
# peak_merge = GenomicRanges::GRanges(readPeakFile(file_bed))


# ##* PEAK ANNOTATION
# ##*
# peakAnno <- annotatePeak(peak_merge, tssRegion=c(-1500, 1500),
#                          TxDb=txdb, annoDb="org.Hs.eg.db")
# p = as.data.frame(peakAnno@anno)

# write.table(p, file="/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/overlap_peaks/FTO_PCDH_overlap_chipseeker_peakAnno.txt", quote = F, row.names = F, sep = "\t")

# plotAnnoBar(peakAnno)

# ##* pathways
# genes = as.data.frame(peakAnno)$geneId
# names(genes) = sub("_", "\n", names(genes))
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            pvalueCutoff  = 0.01,
#                            pAdjustMethod = "BH")
# dotplot(compKEGG, showCategory = 12, title = "KEGG Pathway Enrichment Analysis")

# enrichplot::cnetplot(compKEGG, showCategory = 5, layout = "kk")

# ##********************************************************************************
# ##********************************************************************************

# genelist = list(up = up.1$entrezid, down = down.1$entrezid)
# genelist = (unlist(genelist))
# genelist = sort(genelist, decreasing = FALSE)

# # compKEGG <- compareCluster(geneCluster   = genelist,
# #                            fun           = "enrichPathway",
# #                            pvalueCutoff  = 0.01,
# #                            pAdjustMethod = "BH")

# goa <- limma::goana(genes, species="Hs")

# topGO(goa)

# library(TimeSeriesExperiment)
# TimeSeriesExperiment::plotEnrichment(enrich = goa, n_max = 10)


##----------------------------------------------------------------------##
##-------------------------------- started here --------------------------------------##

# ##* get files in 
files = Sys.glob("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/peaks/genrich/narrowPeak/*.narrowPeak")
print(files)
# 
# peak_merge = tidyChromosomes(peak_merge, keep.X = TRUE, keep.Y = TRUE,
#                              keep.M = FALSE,
#                              keep.nonstandard = FALSE,
#                              genome = NULL)
# 
peak=GenomicRanges::GRangesList(DN1=readPeakFile(files[[1]]),
                                DN2=readPeakFile(files[[2]]),
                                DN3=readPeakFile(files[[3]]),
                                DP1=readPeakFile(files[[4]]),
                                DP2=readPeakFile(files[[5]]),
                                DP3=readPeakFile(files[[6]])
                                )

 
# ##* to do - process with bam files containing only chr 1-22, x, y

DN1<-peak$DN1
DN2=peak$DN2
DN3=peak$DN3
DP1=peak$DP1
DP2=peak$DP2
DP3=peak$DP3

# ##* covplot
DN3@seqinfo@seqnames 
p <- covplot(peak, title = "ATACseq Peaks over Chromosomes")
col <- c(DN1='red', DN2='green', DN3='blue', DP1='pink', DP2='black', DP3='orange')

##* save plot
pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/covplot.pdf")
p
dev.off()


##* profile of atac regions binding to TSS

# tagMatrix <- getTagMatrix(FTO1, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red", title="ATAC binding to TSS regions - FTO1")
# dev.off()
# 
# tagMatrix <- getTagMatrix(FTO2, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="green", title="ATAC binding to TSS regions - FTO2")
# 
# tagMatrix <- getTagMatrix(FTO3, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="pink", title="ATAC binding to TSS regions - FTO3")
# 
# tagMatrix <- getTagMatrix(PCDH1, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="orange", title="ATAC binding to TSS regions - PCDH1")
# 
# tagMatrix <- getTagMatrix(PCDH2, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="blue", title="ATAC binding to TSS regions - PCDH2")
# 
# tagMatrix <- getTagMatrix(PCDH3, windows=promoter)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="purple", title="ATAC binding to TSS regions - PCDH3")

##* avg plots for all cases

promoter <- getPromoters(TxDb=txdb, upstream=500, downstream=500)
tagMatrixList = list(getTagMatrix(peak$DN1, windows=promoter),
                     getTagMatrix(peak$DN2, windows=promoter),
                     getTagMatrix(peak$DN3, windows=promoter),
                     getTagMatrix(peak$DP1, windows=promoter),
                     getTagMatrix(peak$DP2, windows=promoter),
                     getTagMatrix(peak$DP3, windows=promoter))
names(tagMatrixList) = c("DN1", "DN2", "DN3", "DP1", "DP2", "DP3")

##* avg profile of atac peaks to tss region
##* conf is confidence interval
pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/avgprof.pdf")
plotAvgProf(tagMatrixList, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
dev.off()

pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/peak_heatmap.pdf")
tagHeatmap(tagMatrixList, xlim=c(-500, 500), color=NULL)
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------

# pcdh.uni <- read.delim("/Users/pallavi/Documents/cl/projects/atac/FTO_PCDH/data/peaks_diffbind/peaks_PCDH_unique.csv", sep=" ")
# pcdh.uni = pcdh.uni %>% separate(chr.start.end.score, c("chr", "start", "end", "score"))

#---------------------------------------------------------------------------------------------------------------------------------------


##* PEAK ANNOTATION
##*
peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db",
                       verbose=FALSE)

##* write peak files
setwd("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/peaks/genrich/chipseeker")
file_name <- c("DN1.csv", "DN2.csv", "DN3.csv", "DP1.csv", "DP2.csv", "DP3.csv")
lapply(seq_along(peakAnnoList), function(i){
  write.csv(peakAnnoList[[i]]@anno, file_name[i], row.names = FALSE)
})

##* ATAC peak annotation comparison
pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/peak_barplot.pdf")
plotAnnoBar(peakAnnoList)
dev.off()


##* distribution of TF binding loci
pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/peak_dist_to_TSS.pdf")
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci relative to TSS")
dev.off()


##* find peak overlaps - library(ChIPpeakAnno)
peak_overlap_DN <- findOverlapsOfPeaks(peak$DN1,
                                    peak$DN2,
                                    peak$DN3,
                                    maxgap=1000)
saveRDS(peak_overlap_DN, "/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/peaks/genrich/chipseeker/overlap_peaks_DN.rds")
df<-data.frame((peak_overlap_DN$overlappingPeaks))

peak_overlap_DN <- findOverlapsOfPeaks(peak$DP1,
                                    peak$DP2,
                                    peak$DP3,
                                    maxgap=1000)
saveRDS(peak_overlap_DN, "/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/peaks/genrich/chipseeker/overlap_peaks_DP.rds")

##********************************************************************************
##********************************************************************************

##* pathways
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.01,
                           pAdjustMethod = "BH")
pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/pathway_dotplot.pdf")
dotplot(compKEGG, showCategory = 12, title = "KEGG Pathway Enrichment Analysis")
dev.off()

pdf("/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/plots/genrich_chipseeker/pathway_cnetplot.pdf")
enrichplot::cnetplot(compKEGG, showCategory = 5, layout = "kk")
dev.off()


##********************************************************************************
##********************************************************************************

# genelist = list(up = up.1$entrezid, down = down.1$entrezid)
# genelist = (unlist(genelist))
# genelist = sort(genelist, decreasing = FALSE)

# # compKEGG <- compareCluster(geneCluster   = genelist,
# #                            fun           = "enrichPathway",
# #                            pvalueCutoff  = 0.01,
# #                            pAdjustMethod = "BH")

# goa <- limma::goana(genelist, species="Hs")

# topGO(goa)

# library(TimeSeriesExperiment)
# TimeSeriesExperiment::plotEnrichment(enrich = goa, n_max = 10)



#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------


# # peak annotate for unique peaks using chippeak anno

# pcdh.uni <- read.delim("/Users/pallavi/Documents/cl/projects/atac/FTO_PCDH/data/peaks_diffbind/peaks_PCDH_unique.bed", sep="\t")
# pcdh.uni = pcdh.uni %>% na.omit() %>% unique()
# pcdh.uni = unique(pcdh.uni)

# bed2grange = function(i){
#   GRanges(seqnames = i$chr,ranges = IRanges(start = i$start,
#                                             end = i$end,
#                                             names = rownames(i)))
# }

# pcdh.uni.bed = bed2grange(i=pcdh.uni)

# pcdh.uni.peaksAnno=annotatePeak(pcdh.uni.bed, tssRegion=c(-1500, 1500), TxDb=txdb, annoDb="org.Hs.eg.db")
# # Annotated peaks generated by ChIPseeker
# # 11378/11388  peaks were annotated
# # Genomic Annotation Summary:
# #   Feature  Frequency
# # 9           Promoter 23.1675163
# # 4             5' UTR  0.1933556
# # 3             3' UTR  2.2851116
# # 1           1st Exon  1.3007558
# # 7         Other Exon  2.7333451
# # 2         1st Intron 15.1256811
# # 8       Other Intron 27.7025839
# # 6 Downstream (<=300)  0.1757778
# # 5  Distal Intergenic 27.3158727

# write.table(pcdh.uni.peaksAnno, file = paste0(RES, "data/annotated_peaks/Unique_PCDH_annotated.csv"), sep="\t", quote=F, row.names=F)

# ##* Similarly for fto gives-

# # Annotated peaks generated by ChIPseeker
# # 11485/11497  peaks were annotated
# # Genomic Annotation Summary:
# #   Feature  Frequency
# # 9           Promoter 22.3247714
# # 4             5' UTR  0.2002612
# # 3             3' UTR  1.8894210
# # 1           1st Exon  1.2973444
# # 7         Other Exon  2.5163256
# # 2         1st Intron 15.2285590
# # 8       Other Intron 27.7405311
# # 6 Downstream (<=300)  0.1218981
# # 5  Distal Intergenic 28.6808881

# list.annotations <- list(FTO=fto.uni.peaksAnno, PCDH=pcdh.uni.peaksAnno)

# jpeg(file = paste0(RES, "plots/PCDH_FTO_TSS_distribution.jpeg"))
# plotDistToTSS(list.annotations, title="PCDH and FTO - Distribution of transcription factor-binding loci relative to TSS")
# dev.off()

