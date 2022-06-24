
##* packages
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(ggplot2)
library(GenomicRanges)
library(ReactomePA)
library(BRGenomics)


###* path
# RES="/projects/b1017/pallavi/ATAC_Seq/FTO_PCDH/results"
RES = "/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results"
setwd(RES)

##* get files
file_bed = "/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/overlap_peaks/bed_FTO_PCDH_overlap.bed"
peak_merge = GenomicRanges::GRanges(readPeakFile(file_bed))


##* PEAK ANNOTATION
##*
peakAnno <- annotatePeak(peak_merge, tssRegion=c(-1500, 1500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
p = as.data.frame(peakAnno@anno)

write.table(p, file="/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/overlap_peaks/FTO_PCDH_overlap_chipseeker_peakAnno.txt", quote = F, row.names = F, sep = "\t")

plotAnnoBar(peakAnno)

##* pathways
genes = as.data.frame(peakAnno)$geneId
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.01,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 12, title = "KEGG Pathway Enrichment Analysis")

enrichplot::cnetplot(compKEGG, showCategory = 5, layout = "kk")

##********************************************************************************
##********************************************************************************

genelist = list(up = up.1$entrezid, down = down.1$entrezid)
genelist = (unlist(genelist))
genelist = sort(genelist, decreasing = FALSE)

# compKEGG <- compareCluster(geneCluster   = genelist,
#                            fun           = "enrichPathway",
#                            pvalueCutoff  = 0.01,
#                            pAdjustMethod = "BH")

goa <- limma::goana(genes, species="Hs")

topGO(goa)

library(TimeSeriesExperiment)
TimeSeriesExperiment::plotEnrichment(enrich = goa, n_max = 10)


##----------------------------------------------------------------------##
##----------------------------------------------------------------------##

# ##* get files in 
# files = Sys.glob("/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/genrich_peakcalls/*.narrowPeak")
# print(files)
# 
# file_bed = Sys.glob("/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/peaks/*.bed")
# print(file_bed)
# peak_merge = GenomicRanges::GRanges(readPeakFile(file_bed[[1]]))
# peak_merge = tidyChromosomes(peak_merge, keep.X = TRUE, keep.Y = TRUE,
#                              keep.M = FALSE,
#                              keep.nonstandard = FALSE,
#                              genome = NULL)
# 
# peak=GenomicRanges::GRangesList(FTO1=readPeakFile(files[[1]]),
#                                 FTO2=readPeakFile(files[[2]]),
#                                 # FTO3=readPeakFile(files[[3]]),
#                                 PCDH1=readPeakFile(files[[4]]),
#                                 PCDH2=readPeakFile(files[[5]])
#                                 # PCDH3=readPeakFile(files[[6]])
#                                 )
# ##* all the peaks of the exp
# all_peaks = peak
# 
# ##* PC plot shows FTO3 and PCDH3 cluster so remove them and do analysis
# peak$FTO3 <- NULL
# peak$PCDH3 <- NULL
# 
# 
# ##* to do - process with bam files containing only chr 1-22, x, y
# peak = tidyChromosomes(peak, keep.X = TRUE, keep.Y = TRUE,
#                        keep.M = FALSE,
#                        keep.nonstandard = FALSE,
#                        genome = NULL)
# 
# FTO1<-peak$FTO1
# FTO2=peak$FTO2
# # FTO3=peak$FTO3
# PCDH1=peak$PCDH1
# PCDH2=peak$PCDH2
# # PCDH3=peak$PCDH3
# 
# ##* covplot
# # FTO1@seqinfo@seqnames 
# p <- covplot(peak, title = "ATACseq Peaks over Chromosomes")
# # col <- c(FTO1='red', FTO2='green', FTO3='blue', PCDH1='pink', PCDH2='black', PCDH3='orange')
# col <- c(FTO1='red', FTO2='green', PCDH1='pink', PCDH2='black')
# 
# ##* save plot
# pdf("/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/plots/covplot.pdf") 
# p
# # p + facet_grid(chr ~ .id) + scale_color_manual(values=col) + scale_fill_manual(values=col)
# dev.off() 
# 
# ##* profile of atac regions binding to TSS
# # promoter1 <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
# 
# tagMatrix <- getTagMatrix(FTO1, windows=promoter)
# # pdf("/Users/pallavi/Documents/q/ATAC_Seq/FTO_PCDH/results/plots/tagheatmap_fto1.pdf",
# #     width = 12, height = 7) 
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
# 
# ##* avg plots for all cases
# ##* 
# promoter <- getPromoters(TxDb=txdb, upstream=500, downstream=500)
# tagMatrixList = list(getTagMatrix(peak$FTO1, windows=promoter),
#                      getTagMatrix(peak$FTO2, windows=promoter),
#                      getTagMatrix(peak$PCDH1, windows=promoter),
#                      getTagMatrix(peak$PCDH2, windows=promoter))
# 
# 
# ##* avg profile of atac peaks to tss region
# ##* conf is confidence interval
# plotAvgProf(tagMatrixList, xlim=c(-500, 500),
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf=0.95)
# 
# plotAvgProf(tagMatrixList, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
# 
# # tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
# 
# 
# 
# ##* PEAK ANNOTATION
# ##* 
# peakAnno <- annotatePeak(peak_merge, tssRegion=c(-3000, 3000),
#                          TxDb=txdb, annoDb="org.Hs.eg.db")
# write.table(as.data.frame(peakAnno@anno), file="merged_chipseeker_peakAnno.txt", quote = F, row.names = F, sep = "\t")
# 
# p = as.data.frame(peakAnno@anno)
# 
# 
# peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb,
#                        tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db",
#                        verbose=FALSE)
# 
# peakAnnoList1 = peakAnnoList
# 
# ##* write peak files
# file_name <- c("FTO1.csv", "FTO2.csv", "FTO3.csv", "PCDH1.csv", "PCDH2.csv", "PCDH3.csv")
# lapply(seq_along(peakAnnoList), function(i){
#   write.csv(peakAnnoList[[i]]@anno, file_name[i], row.names = FALSE)
# })
# 
# ##* ATAC peak annotation comparison
# plotAnnoBar(peakAnnoList)
# 
# 
# ##* distribution of TF binding loci
# plotDistToTSS(peakAnnoList,
#               title="Distribution of transcription factor-binding loci relative to TSS")
# 
# 
# ##* find peak overlaps - library(ChIPpeakAnno)
# peak_overlap <- findOverlapsOfPeaks(peak$FTO1, 
#                                     peak$FTO2, 
#                                     peak$PCDH1,
#                                     peak$PCDH2,
#                                     maxgap=1000)
# saveRDS(peak_overlap, "overlap_peaks_chipseeker.rds")
# 
# ##********************************************************************************
# ##********************************************************************************
# 
# ##* pathways
# genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) = sub("_", "\n", names(genes))
# compKEGG <- compareCluster(geneCluster   = genes,
#                            fun           = "enrichKEGG",
#                            pvalueCutoff  = 0.01,
#                            pAdjustMethod = "BH")
# dotplot(compKEGG, showCategory = 12, title = "KEGG Pathway Enrichment Analysis")
# 
# enrichplot::cnetplot(compKEGG, showCategory = 5, layout = "kk")
# 
# ##********************************************************************************
# ##********************************************************************************
# 
# genelist = list(up = up.1$entrezid, down = down.1$entrezid)
# genelist = (unlist(genelist))
# genelist = sort(genelist, decreasing = FALSE)
# 
# # compKEGG <- compareCluster(geneCluster   = genelist,
# #                            fun           = "enrichPathway",
# #                            pvalueCutoff  = 0.01,
# #                            pAdjustMethod = "BH")
# 
# goa <- limma::goana(genelist, species="Hs")
# 
# topGO(goa)
# 
# library(TimeSeriesExperiment)
# TimeSeriesExperiment::plotEnrichment(enrich = goa, n_max = 10)
# 
# 
# 
