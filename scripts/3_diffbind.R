#---------------------------------------------------------------------------------------------------------------------------------------

# https://www.one-tab.com/page/JUc6Sh1MTFuqRe6uyyi74A

#---------------------------------------------------------------------------------------------------------------------------------------
#- conda activate R411

##* packages needed
library(DiffBind) #DiffBind 3.4 used
library(parallel)
library(dplyr)
library(edgeR)
# library(profileplyr) #profile plots
library(ChIPpeakAnno)
# Loading TSS Annotation For Human Sapiens (GRCh38) Obtained From BiomaRt -- inside Chippeakanno
data(TSS.human.GRCh38)

##* path
RES="/home/psurana/projects/Matei_lab/ofkar_DP_DN/data/res/"

#---------------------------------------------------------------------------------------------------------------------------------------

##* read samplesheet
samplesheet = read.csv(file = paste0(RES, "try.csv"))
samplesheet


##* Read peaksets 
dbObj <- dba(sampleSheet=samplesheet)
dbObj
# ExptConsensus <-  dba(dbObj, mask=dbObj$masks$Consensus)
# ConsensusPeaks <- dba.peakset(dbObj, bRetrieve=TRUE)
#- macs2
# 6 Samples, 106611 sites in matrix (296003 total):
#   ID   Tissue Factor Condition Replicate Intervals
# 1  1 multiple     DN  stemness         1    117640
# 2  2 multiple     DN  stemness         2    160460
# 3  3 multiple     DN  stemness         3    142267
# 4  4 multiple     DP   unknown         1    133371
# 5  5 multiple     DP   unknown         2    182011
# 6  6 multiple     DP   unknown         3     67915
#- genrich
# 6 Samples, 73347 sites in matrix (112901 total):
#   ID   Tissue Factor Condition Replicate Intervals
# 1  1 multiple     DN        DN         1     69179
# 2  2 multiple     DN        DN         2     75508
# 3  3 multiple     DN        DN         3     69892
# 4  4 multiple     DP        DP         1     57717
# 5  5 multiple     DP        DP         2     70892
# 6  6 multiple     DP        DP         3     65263
#- genrich - p=0.1
# 6 Samples, 133893 sites in matrix (248338 total):
#   ID   Tissue Factor Condition Replicate Intervals
# 1  1 multiple     DN        DN         1    129005
# 2  2 multiple     DN        DN         2    127602
# 3  3 multiple     DN        DN         3    151469
# 4  4 multiple     DP        DP         1    112962
# 5  5 multiple     DP        DP         2    124168
# 6  6 multiple     DP        DP         3    152279

##* affinity binding matrix - - count read mapping to the peaks
##- summits=100 for ATAC-seq, which results in 201bp windows
#-- READS- (the "Full" library sizes). 
#-- FRiP (Fraction of Reads in Peaks) - This is the proportion of reads for that sample that overlap a peak in the consensus peakset, and can be used to indicate which samples show more enrichment overall.
#-- For each sample, multiplying the value in the Reads column by the corresponding FRiP value will yield the number of reads that overlap a consensus peak. T

dbObj <- dba.count(dbObj, summits=50)
dbObj
# 6 Samples, 102276 sites in matrix: #- macs2
#   ID   Tissue Factor Condition Replicate    Reads FRiP
# 1  1 multiple     DN  stemness         1 34995189 0.06
# 2  2 multiple     DN  stemness         2 38734649 0.08
# 3  3 multiple     DN  stemness         3 28693634 0.06
# 4  4 multiple     DP   unknown         1 39506275 0.05
# 5  5 multiple     DP   unknown         2 42316137 0.07
# 6  6 multiple     DP   unknown         3 23651148 0.05
# 6 Samples, 70927 sites in matrix: #- genrich
#   ID   Tissue Factor Condition Replicate    Reads FRiP
# 1  1 multiple     DN        DN         1 34995189 0.04
# 2  2 multiple     DN        DN         2 38734649 0.06
# 3  3 multiple     DN        DN         3 28693634 0.04
# 4  4 multiple     DP        DP         1 39506275 0.04
# 5  5 multiple     DP        DP         2 42316137 0.05
# 6  6 multiple     DP        DP         3 23651148 0.04
# 6 Samples, 128458 sites in matrix: #- genrich = p = 0.1
#   ID   Tissue Factor Condition Replicate    Reads FRiP
# 1  1 multiple     DN        DN         1 34995189 0.05
# 2  2 multiple     DN        DN         2 38734649 0.07
# 3  3 multiple     DN        DN         3 28693634 0.05
# 4  4 multiple     DP        DP         1 39506275 0.05
# 5  5 multiple     DP        DP         2 42316137 0.06
# 6  6 multiple     DP        DP         3 23651148 0.05
dbObj <- dba.count(dbObj, summits=0) #- genrich = p = 0.1 (summits=0)
dbObj
# 4 Samples, 98673 sites in matrix:
#   ID   Tissue Factor Condition Replicate    Reads FRiP
# 1  1 multiple     DN        DN         1 34995189 0.26
# 2  2 multiple     DN        DN         2 38734649 0.29
# 3  3 multiple     DP        DP         1 39506275 0.23
# 4  4 multiple     DP        DP         2 42316137 0.28
#---------------------------------------------------------------------------------------------------------------------------------------

##* plots
jpeg(file = paste0(RES, "plots/genirch_p=0.1_summits=0_diffbind.jpeg"))
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()


#correlation between samples plot
jpeg(file = paste0(RES, "plots/genirch_p=0.1__summits=0corr_diffbind.jpeg"))
dba.plotHeatmap(dbObj)
dev.off()

##* normalize lib sizes for fair comparision across diff library sizes
dbObj <- dba.normalize(dbObj)
dbObj

#---------------------------------------------------------------------------------------------------------------------------------------

##* which samples to compare?
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)
dbObj
# Design: [~Factor] | 1 Contrast:
#   Factor Group Samples Group2 Samples2
# 1 Factor    DP       3     DN        3


##* differential enrichment analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS, bBlacklist = FALSE, bGreylist = FALSE)
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)

##* inspect obj
dba.show(dbObj, bContrasts=T)
#   Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2 #- peak calling with macs2
# 1 Factor    DP       3     DN        3        0         0

#---------------------------------------------------------------------------------------------------------------------------------------

##* VIZ

#plot corr heatmap of 1320 differential sites
jpeg(file = paste0(RES, "plots/heatmap_diff_sites_0_Diffbind.jpeg"))
plot(dbObj, contrast=1)
dev.off()


#significant regions dea pc plot
jpeg(file = paste0(RES, "plots/pca_diff_anal_Diffbind.jpeg"))
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

# jpeg(file = paste0(RES, "plots/venn_diffbind.jpeg"))
# dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
# dev.off()

#overlap between peaks
jpeg(file = paste0(RES, "plots/peak_overlap_venn_diff_anal_Diffbind.jpeg"))
dba.plotVenn(dbObj, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()

jpeg(file = paste0(RES, "plots/maPlot_diffbind.jpeg"))
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()

jpeg(file = paste0(RES, "plots/BoxPlot_diffbind.jpeg"))
dba.plotBox(dbObj)
dev.off()

pvals = dba.plotBox(dbObj)


#volcano plot
jpeg(file = paste0(RES, "plots/volcanoplot_diffbind.jpeg"))
dba.plotVolcano(dbObj)
dev.off()

#profile plot - merge all samples of a type as a condition
profiles <- dba.plotProfile(dbObj,merge=c(DBA_TISSUE, DBA_REPLICATE))
dba.plotProfile(profiles)

#---------------------------------------------------------------------------------------------------------------------------------------

##* results
results <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# sum(results$Fold>0)
#[1] 38863
# 26800 -- genrich
# sum(results$Fold<0)
#[1] 62482
# 43623 -- genrich

##* binding site overlap plot pie -- not needed -- only one contrast
# jpeg(file = paste0(RES, "plots/binding_site_overlap_diffbind.jpeg"))
# dba.plotVenn(dbObj, contrast=1)
# dev.off()

out <- as.data.frame(results)
write.table(out, file = paste0(RES, "diffbind/diffbind_samples.csv"), sep="\t", quote=F, row.names=F)

#---------------------------------------------------------------------------------------------------------------------------------------

##* save session
save.image(paste0(RES, "data/diffBind_FTO_PCDH.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------

##* differentially expressed sites -- not working - fdr cutoff not working here
dp <- out %>% 
  filter(Fold > 0)
write.table(dp, file = paste0(RES, "diffbind/dp_differential_sites.csv"), sep="\t", quote=F, row.names=F)

dn <- out %>% 
  filter(Fold < 0)
#[1] 76 11
write.table(dn, file = paste0(RES, "diffbind/dn_differential_sites.csv"), sep="\t", quote=F, row.names=F)

#---------------------------------------------------------------------------------------------------------------------------------------

##* overlap analysis
olap.rate <- dba.overlap(dbObj,mode=DBA_OLAP_RATE)
olap.rate 
#[1] 100342 100342  92070  84888  77581  67598



#---------------------------------------------------------------------------------------------------------------------------------------

##* save session
save.image(paste0(RES, "diffbind/diffBind_DP_DN.RData"))



