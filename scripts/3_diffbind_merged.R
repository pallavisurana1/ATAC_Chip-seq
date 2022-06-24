


#---------------------------------------------------------------------------------------------------------------------------------------

# https://www.one-tab.com/page/JUc6Sh1MTFuqRe6uyyi74A
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/diffBind/lab-diffBinding-remote.html
# nohup R CMD BATCH 3_diffbind.R &> out/7_diffbind.out &

#---------------------------------------------------------------------------------------------------------------------------------------
#- conda activate R411


##* packages needed
pacman::p_load(parallel, dplyr, BiocManager, ReactomePA, clusterProfiler, biomaRt)
pacman::p_load(DiffBind, edgeR, ChIPpeakAnno,ChIPseeker, ggplot2)
data(TSS.human.GRCh38)

##* path
RES="/home/psurana/projects/Matei_lab/chip_cut_run/data/res/"
# RES="Z:/projects/Matei_lab/chip_cut_run/data/res/"
# RES="/Users/pallavi/Documents/m/projects/Matei_lab/chip_cut_run/data/res/"

#---------------------------------------------------------------------------------------------------------------------------------------

##* read samplesheet
samplesheet = read.csv(file = paste0(RES, "diffbind_merged/4_samples.csv"), row.names=F)
samplesheet %>% head()


##* Read peaksets 
dbObj <- dba(sampleSheet=samplesheet)
dbObj
## ExptConsensus <-  dba(dbObj, mask=dbObj$masks$Consensus)
## ConsensusPeaks <- dba.peakset(res.dbObj, bRetrieve=TRUE)


##* affinity binding matrix - - count read mapping to the peaks
##- summits=100 for ATAC-seq, which results in 201bp windows
#-- READS- (the "Full" library sizes). 
#-- FRiP (Fraction of Reads in Peaks) - This is the proportion of reads for that sample that overlap a peak in the consensus peakset, and can be used to indicate which samples show more enrichment overall.
#-- For each sample, multiplying the value in the Reads column by the corresponding FRiP value will yield the number of reads that overlap a consensus peak. T

# dbObj <- dba.count(dbObj, summits=50, config=data.frame(RunParallel=TRUE))
dbObj = dba.count(dbObj, minOverlap=3, score=DBA_SCORE_TMM_MINUS_FULL, fragmentSize=130)
dbObj



##* normalize lib sizes for fair comparision across diff library sizes -TMM
dbObj=dba.normalize(dbObj, normalize=DBA_NORM_TMM)
dbObj
#- inspect normalization
dba.normalize(dbObj, bRetrieve=TRUE)

# dbObj$class[DBA_TISSUE, ] = rep(c("DN", "DP"), each =12)


##* which samples to compare?
dbObj = dba.contrast(dbObj, categories=c(DBA_TISSUE, DBA_FACTOR), minMembers=4)

# inspecting the object: how many contrasts were set in the previous step
dbObj
#---------------------------------------------------------------------------------------------------------------------------------------

#- some plots of data exploration
#- plotting the correlation of libraries based on normalised counts of reads in peaks
jpeg(file = paste0(RES, "diffbind_merged/plots/correlation_libraries_normalised.png"))
plot(dbObj)
dev.off()

# PCA scores plot: data overview
pdf(file = paste0(RES, "diffbind_merged/plots/macs2_PCA_diffbind_all.pdf"),  width = 9, height = 7)
plot = dba.plotPCA(dbObj,DBA_FACTOR,label=DBA_CONDITION, dotSize=1, labelSize=0.4)
print(plot)
dev.off()


#--------------------- DEA ------------------------------------------------------------------------------------------------------------------
load(paste0(RES, "diffbind_merged/diffBind_DP_DN.RData"))

##* differential enrichment analysis
res.dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)

save.image(paste0(RES, "diffbind_merged/diffBind_DP_DN.RData"))
# load(paste0(RES, "diffbind_merged/diffBind_DP_DN.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------

# #* VIZ

# #plot corr heatmap of 1320 differential sites
# pdf(file = paste0(RES, "diffbind_merged/plots/heatmap_diff_sites_Diffbind.pdf"))
# plot(res.dbObj, contrast=3)
# plot(res.dbObj, contrast=8)
# plot(res.dbObj, contrast=12)
# dev.off()

# #significant regions dea pc plot
# pdf(file = paste0(RES, "diffbind_merged/plots/pca_diff_anal_Diffbind.pdf"))
# dba.plotPCA(res.dbObj, contrast=3, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
# dba.plotPCA(res.dbObj, contrast=8, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
# dba.plotPCA(res.dbObj, contrast=12, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
# dev.off()



# # jpeg(file = paste0(RES, "plots/maPlot_diffbind.jpeg"))
# # dba.plotMA(res.dbObj, contrast = c(3,8,12), method=DBA_DESEQ2)
# # dev.off()

# # jpeg(file = paste0(RES, "plots/BoxPlot_diffbind.jpeg"))
# # dba.plotBox(res.dbObj, contrast = c(3,8,12), method=DBA_DESEQ2)
# # dev.off()

# # pvals = dba.plotBox(dbObj)






# #---------------------------------------------------------------------------------------------------------------------------------------

# #---------------------------------------------------------------------------------------------------------------------------------------

##- diff expr for consensus peaks
dbObj.peaks <- dba.count(res.dbObj, peaks=res.dbObj$masks$Consensus, summits=T)

##* normalize lib sizes for fair comparision across diff library sizes
dbObj.peaks <- dba.normalize(dbObj.peaks)
dbObj.peaks

##* which samples to compare?
dbObj.peaks <- dba.contrast(dbObj.peaks, categories=DBA_FACTOR, minMembers = 4)
dbObj.peaks


##* differential enrichment analysis
dbObj.peaks <- dba.analyze(dbObj.peaks, method=DBA_DESEQ2)

##* inspect obj
dba.show(dbObj.peaks, bContrasts = T)

# #volcano plot
pdf(file = paste0(RES, "diffbind_merged/plots/volcanoplot_diffbind.pdf"))
dba.plotVolcano(dbObj.peaks, contrast = 3, method=DBA_DESEQ2)
dba.plotVolcano(dbObj.peaks, contrast = 8, method=DBA_DESEQ2)
dba.plotVolcano(dbObj.peaks, contrast = 12, method=DBA_DESEQ2)
dev.off()

#profile plot - merge all samples of a type as a condition
profiles <- dba.plotProfile(dbObj.peaks, samples=list(DP=dbObj.peaks$mask$DP,
                                                      DN=dbObj.peaks$mask$DN))
pdf(file = paste0(RES, "diffbind_merged/plots/profileplot_diffbind.pdf"), width=9, height=9)
dba.plotProfile(profiles)
dev.off()

#DNK27AC
profiles <- dba.plotProfile(dbObj.peaks, samples=list(DN=dbObj.peaks$mask$DN_K27AC,
                                                      DP=dbObj.peaks$mask$DP_K27AC))
pdf(file = paste0(RES, "diffbind_merged/plots/profileplot_diffbind_DN_DP_K27ac.pdf"), width=9, height=9)
dba.plotProfile(profiles, top_anno_height = unit(1, "cm"))
dev.off()


profiles <- dba.plotProfile(dbObj.peaks, samples=list(DN=dbObj.peaks$mask$DN_K27ME3,
                                                      DP=dbObj.peaks$mask$DP_K27ME3))
pdf(file = paste0(RES, "diffbind_merged/plots/profileplot_diffbind_DN_DP_K27ME3.pdf"), width=9, height=9)
dba.plotProfile(profiles, top_anno_height = unit(1, "cm"))
dev.off()


profiles <- dba.plotProfile(dbObj.peaks, samples=list(DN=dbObj.peaks$mask$DN_K4ME3,
                                                      DP=dbObj.peaks$mask$DP_K4ME3))
pdf(file = paste0(RES, "diffbind_merged/plots/profileplot_diffbind_DN_DP_K4ME3.pdf"), width=9, height=9)
dba.plotProfile(profiles, top_anno_height = unit(1, "cm"))
dev.off()


##overlap between peaks
pdf(file = paste0(RES, "diffbind_merged/plots/peak_overlap_venn_diff_anal_Diffbind.pdf"))
dba.plotVenn(dbObj.peaks, contrast = c(3,8,12), method=DBA_DESEQ2)
dev.off()
#---------------------------------------------------------------------------------------------------------------------------------------


##* differentially expressed sites
#- contrast 3,8,12
results.peaks <- dba.report(dbObj.peaks, method=DBA_DESEQ2, contrast = 8, th=1)
out = as.data.frame(results.peaks)

write.table(out, file = paste0(RES, "diffbind_merged/consensus_peaks/diffbind_res_k27ac.csv"), sep="\t", quote=F, row.names=F)
write.table(out, file = paste0(RES, "diffbind_merged/consensus_peaks/diffbind_res_k27me3.csv"), sep="\t", quote=F, row.names=F)
write.table(out, file = paste0(RES, "diffbind_merged/consensus_peaks/diffbind_res_k4me3.csv"), sep="\t", quote=F, row.names=F)

sum(results.peaks$Fold>0 & results.peaks$FDR<0.05 & results.peaks$'p-value'<0.01)
sum(results.peaks$Fold<0 & results.peaks$FDR<0.05 & results.peaks$'p-value'<0.01)


x <- out %>% filter(FDR<0.05 & Fold>0 & p.value<0.01)
write.table(x, file = paste0(RES, "diffbind_merged/consensus_peaks/fdr_0.05/0.05_k27me3_differential_sites_gain.csv"), sep="\t", quote=F, row.names=F)

x <- out %>% filter(FDR<0.01 & Fold>0 & p.value<0.01)
write.table(x, file = paste0(RES, "diffbind_merged/consensus_peaks/fdr_0.01/0.01_k27me3_differential_sites_gain.csv"), sep="\t", quote=F, row.names=F)

y <- out %>% filter(FDR<0.05 & Fold<0 & p.value<0.01)
write.table(y, file = paste0(RES, "diffbind_merged/consensus_peaks/fdr_0.05/0.05_k27me3_differential_sites_loss.csv"), sep="\t", quote=F, row.names=F)

y <- out %>% filter(FDR<0.01 & Fold<0 & p.value<0.01)
write.table(y, file = paste0(RES, "diffbind_merged/consensus_peaks/fdr_0.01/0.01_k27me3_differential_sites_loss.csv"), sep="\t", quote=F, row.names=F)

dim(x)
dim(y)
#---------------------------------------------------------------------------------------------------------------------------------------

##* save session
save.image(paste0(RES, "diffbind_merged/diffBind_DP_DN.RData"))
saveRDS(res.dbObj, paste0(RES, "diffbind_merged/diffBind_res_only.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------




# > dba.show(dbObj.peaks)
# ID Tissue    Factor    Condition Replicate Caller Intervals    Reads FRiP
# 1   1     DN  DN_K27AC  DN1_K27AC_1         1 counts   2053177  9968048 0.74
# 2   2     DN  DN_K27AC  DN1_K27AC_2         2 counts   2053177  9850672 0.74
# 3   3     DN DN_K27ME3 DN1_K27ME3_1         1 counts   2053177 12103192 0.74
# 4   4     DN DN_K27ME3 DN1_K27ME3_2         2 counts   2053177 11974843 0.74
# 5   5     DN  DN_K4ME3  DN1_K4ME3_1         1 counts   2053177 12601229 0.72
# 6   6     DN  DN_K4ME3  DN1_K4ME3_2         2 counts   2053177 12455977 0.72
# 7   7     DN  DN_K27AC  DN2_K27AC_1         3 counts   2053177 10449450 0.68
# 8   8     DN  DN_K27AC  DN2_K27AC_2         4 counts   2053177 10259731 0.68
# 9   9     DN DN_K27ME3 DN2_K27ME3_1         3 counts   2053177  7249690 0.76
# 10 10     DN DN_K27ME3 DN2_K27ME3_2         4 counts   2053177  7120523 0.76
# 11 11     DN  DN_K4ME3  DN2_K4ME3_1         3 counts   2053177 10404286 0.74
# 12 12     DN  DN_K4ME3  DN2_K4ME3_2         4 counts   2053177 10224253 0.74
# 13 13     DP  DP_K27AC  DP1_K27AC_1         1 counts   2053177  5406523 0.75
# 14 14     DP  DP_K27AC  DP1_K27AC_2         2 counts   2053177  5401076 0.75
# 15 15     DP DP_K27ME3 DP1_K27ME3_1         1 counts   2053177  4451881 0.78
# 16 16     DP DP_K27ME3 DP1_K27ME3_2         2 counts   2053177  4445187 0.79
# 17 17     DP  DP_K4ME3  DP1_K4ME3_1         1 counts   2053177 10219801 0.73
# 18 18     DP  DP_K4ME3  DP1_K4ME3_2         2 counts   2053177 10309458 0.73
# 19 19     DP  DP_K27AC  DP2_K27AC_1         3 counts   2053177 10017957 0.74
# 20 20     DP  DP_K27AC  DP2_K27AC_2         4 counts   2053177  9901244 0.74
# 21 21     DP DP_K27ME3 DP2_K27ME3_1         3 counts   2053177  4724146 0.76
# 22 22     DP DP_K27ME3 DP2_K27ME3_2         4 counts   2053177  4672608 0.76
# 23 23     DP  DP_K4ME3  DP2_K4ME3_1         3 counts   2053177 14261602 0.72
# 24 24     DP  DP_K4ME3  DP2_K4ME3_2         4 counts   2053177 14087661 0.72



