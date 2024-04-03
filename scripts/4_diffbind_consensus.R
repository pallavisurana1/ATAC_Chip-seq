

##* What is FRIP?
# The second is labeled FRiP, which stands for Fraction of Reads in Peaks. 
# This is the proportion of reads for that sample that overlap a peak in the consensus peakset, and can be used to indicate which samples show more enrichment overall.
#---------------------------------------------------------------------------------------------------------------------------------------

##* packages needed
library(DiffBind)
library(parallel)
library(dplyr)
library(edgeR)
# library(profileplyr) #profile plots
library(ChIPpeakAnno)
# Loading TSS Annotation For Human Sapiens (GRCh38) Obtained From BiomaRt -- inside Chippeakanno
data(TSS.human.GRCh38)

#---------------------------------------------------------------------------------------------------------------------------------------

##* path
RES="/home/psurana/projects/atac/FTO_PCDH/"

##* load diffbind image
load(paste0(RES, "data/diffBind_FTO_PCDH.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------

##* Choosing the peaks for the interesting comparison

#get contrasts
dba.show(dbObj, bContrast=T)

data.peaks = dba.report(dbObj, contrast=1)
head(data.peaks)

##* save peaks obj
save(data.peaks, file = paste0(RES, "data/diffBind_peaks_FTO_PCDH.RData"))
load(file = paste0(RES, "data/diffBind_peaks_FTO_PCDH.RData"))



##* gain of func and loss of func
# data.peaks.gain = dba.report(dbObj, bGain=TRUE)
# data.peaks.loss = dba.report(dbObj, bLoss=TRUE)


##* annotate peaks
data.peaksAnno=annotatePeakInBatch(data.peaks, AnnotationData=TSS.human.GRCh38)
#1320 x 15 metadata cols
write.table(data.peaksAnno, file = paste0(RES, "data/annotated_peaks/peaks_FTO_vs_PCDH_annotated.csv"), sep="\t", quote=F, row.names=F)

##* annotate all results
results.peaksAnno=annotatePeakInBatch(results, AnnotationData=TSS.human.GRCh38)
#100292 x 15 metadata cols
write.table(results.peaksAnno, file = paste0(RES, "data/annotated_peaks/peaks_all_annotated.csv"), sep="\t", quote=F, row.names=F)

#-------------- or ------------------


##* get consensus peaks
db_consensus <- dba(dbObj, mask=dbObj$masks$Consensus, minOverlap=2)

# FTO_consensus = dba.overlap(db_consensus, c(1,2,3),mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
# PCDH_consensus = dba.overlap(db_consensus, c(4,5,6),mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)

db_overlaps  <- dba.overlap(db_consensus, mask=1:2, DataType=DBA_DATA_FRAME)


db_overlaps$inAll %>% dim()   
#[1] 72490     5
db_overlaps$onlyA %>% dim()
#[1] 11509     4
db_overlaps$onlyB %>% dim()
#[1] 11396     4

write.table(db_overlaps$inAll, file = paste0(RES, "data/peaks_diffbind/peaks_FTO_vs_PCDH_olap.csv"), sep="\t", quote=F, row.names=F)
write.table(db_overlaps$onlyA, file = paste0(RES, "data/peaks_diffbind/peaks_FTO_unique.csv"), sep="\t", quote=F, row.names=F)
write.table(db_overlaps$onlyB, file = paste0(RES, "data/peaks_diffbind/peaks_PCDH_unique.csv"), sep="\t", quote=F, row.names=F)


##* annotate consensus peaks
consensus.peaksAnno=annotatePeakInBatch(consensus_peaks, AnnotationData=TSS.human.GRCh38)
#100292 x 15 metadata cols
write.table(consensus.peaksAnno, file = paste0(RES, "data/peaks_FTO_vs_PCDH_annotated.csv"), sep="\t", quote=F, row.names=F)

#---------------------------------------------------------------------------------------------------------------------------------------


##* rerun analysis using consensus peaks
dbObj.peaks <- dba.count(dbObj, peaks=dbObj$masks$Consensus, summits=T)


##* normalize lib sizes for fair comparision across diff library sizes
dbObj.peaks <- dba.normalize(dbObj.peaks)
dbObj.peaks
# 6 Samples, 100224 sites in matrix:
#   ID   Tissue Factor Condition Replicate    Reads FRiP
# 1  1 multiple    FTO  knockout         1 53826318 0.29
# 2  2 multiple    FTO  knockout         2 71188301 0.30
# 3  3 multiple    FTO  knockout         3 76174504 0.31
# 4  4 multiple   PCDH   control         1 80914933 0.35
# 5  5 multiple   PCDH   control         2 55075748 0.32
# 6  6 multiple   PCDH   control         3 70921387 0.32

##* which samples to compare?
dbObj.peaks <- dba.contrast(dbObj.peaks, categories=DBA_FACTOR, minMembers = 2)
dbObj.peaks
# Design: [~Factor] | 1 Contrast:
#   Factor Group Samples Group2 Samples2
# 1 Factor  PCDH       3    FTO       3


##* differential enrichment analysis
dbObj.peaks <- dba.analyze(dbObj.peaks, method=DBA_ALL_METHODS)


##* inspect obj
dba.show(dbObj.peaks, bContrasts = T)
#   Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Factor  PCDH       3    FTO        3        4      1437

#---------------------------------------------------------------------------------------------------------------------------------------


##* differentially expressed sites

results.peaks <- dba.report(dbObj.peaks, method=DBA_DESEQ2, contrast = 1, th=1)
sum(results.peaks$Fold>0 & results.peaks$FDR<0.05 & results.peaks$'p-value'<0.01)
sum(results.peaks$Fold<0 & results.peaks$FDR<0.05 & results.peaks$'p-value'<0.01)

out <- as.data.frame(results.peaks)
# [1] 100224     11
write.table(out, file = paste0(RES, "data/consensus_peaks_diffbind/diffbind_samples.csv"), sep="\t", quote=F, row.names=F)

pcdh <- out %>% filter(FDR<0.05 & Fold>0 & p.value<0.01)
#[1] 1336   11
write.table(pcdh, file = paste0(RES, "data/consensus_peaks_diffbind/PCDH_differential_sites.csv"), sep="\t", quote=F, row.names=F)

fto <- out %>% filter(FDR<0.05 & Fold<0 & p.value<0.01)
#[1] 101 11
write.table(fto, file = paste0(RES, "data/consensus_peaks_diffbind/FTO_differential_sites.csv"), sep="\t", quote=F, row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------

##* save session
save.image(paste0(RES, "data/diffBind_overlaps_FTO_PCDH.RData"))

load(paste0(RES, "data/diffBind_overlaps_FTO_PCDH.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------



