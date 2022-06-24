

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
RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/peaks/macs2/narrowPeak/"
setwd(RES)

## to check peak profiles of individual peaks
files = Sys.glob("*.narrowPeak")
print(files)

peak = lapply(files, readPeakFile)

names_peak = c("DN1_K27AC_1",	"DN1_K27AC_2",	"DN1_K27ME3_1",	
               "DN1_K27ME3_2",	"DN1_K4ME3_1",	"DN1_K4ME3_2",	
               "DN2_K27AC_1",	"DN2_K27AC_2",	"DN2_K27ME3_1",	
               "DN2_K27ME3_2",	"DN2_K4ME3_1",	"DN2_K4ME3_2",	
               "DP1_K27AC_1",	"DP1_K27AC_2",	"DP1_K27ME3_1",	
               "DP1_K27ME3_2",	"DP1_K4ME3_1",	"DP1_K4ME3_2",	
               "DP2_K27AC_1", "DP2_K27AC_2",	"DP2_K27ME3_1",	
               "DP2_K27ME3_2",	"DP2_K4ME3_1",	"DP2_K4ME3_2")

names(peak) = names_peak



RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/"


peak_dn1 = peak[1:6]
peak_dn2 = peak[7:12]
peak_dp1 = peak[13:18]
peak_dp2 = peak[19:24]

#- took too long
# pdf(file = paste0(RES, "diffbind_merged/plots/chipseeker_individual/coverage_plot.pdf"))
# s <- covplot(peak_dn1, 
#              title = "Chipseq Peaks - DN1 over Chromosomes")
# s
# s <- covplot(peak_dn2, 
#              title = "Chipseq Peaks - DN2 over Chromosomes")
# s
# s <- covplot(peak_dp1, 
#              title = "Chipseq Peaks - DP1 over Chromosomes")
# s
# s <- covplot(peak_dp2, 
#              title = "Chipseq Peaks - DP2 over Chromosomes")
# s
# dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=500, downstream=500)

tagmat1 = lapply(peak_dn1, getTagMatrix, windows=promoter)
tagmat2 = lapply(peak_dn2, getTagMatrix, windows=promoter)
tagmat3 = lapply(peak_dp1, getTagMatrix, windows=promoter)
tagmat4 = lapply(peak_dp2, getTagMatrix, windows=promoter)

pdf(file = paste0(RES, "diffbind_merged/plots/chipseeker_individual/avgprofile.pdf"))
plotAvgProf(tagmat1, xlim=c(-500, 500),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagmat1, xlim=c(-500, 500),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagmat2, xlim=c(-500, 500),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagmat3, xlim=c(-500, 500),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagmat4, xlim=c(-500, 500),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagmat1, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
plotAvgProf(tagmat2, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
plotAvgProf(tagmat3, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
plotAvgProf(tagmat4, xlim=c(-500, 500), conf=0.95,resample=500, facet="row")
dev.off()



#- peaks heatmap
pdf(file = paste0(RES, "diffbind_merged/plots/chipseeker_individual/peakheatmap.pdf"))
tagHeatmap(tagmat1, xlim=c(-500, 500), color=rainbow(length(tagmat1)))
tagHeatmap(tagmat2, xlim=c(-500, 500), color=rainbow(length(tagmat2)))
tagHeatmap(tagmat3, xlim=c(-500, 500), color=rainbow(length(tagmat3)))
tagHeatmap(tagmat4, xlim=c(-500, 500), color=rainbow(length(tagmat4)))
dev.off()


