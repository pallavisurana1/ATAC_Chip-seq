
#---------------------------------------------------------------------------------------------------------------------------------------

RES="/Users/pallavi/Documents/c/projects/atac/FTO_PCDH/"
load(file = paste0(RES, "data/genes_list.RData"))

#---------------------------------------------------------------------------------------------------------------------------------------

compKEGG <- compareCluster(geneCluster = genes,
                           fun = "enrichKEGG",
                           organism = "hsa",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 12, title = "KEGG Pathway Enrichment Analysis")

#---------------------------------------------------------------------------------------------------------------------------------------

genelist = (unlist(genes))
genelist = sort(genelist, decreasing = FALSE)

goa.fto <- limma::goana(genes$FTO, species="Hs")
goa.pcdh <- limma::goana(genes$PCDH, species="Hs")

limma::topGO(goa.pcdh)
limma::topGO(goa.fto)

library(TimeSeriesExperiment)
TimeSeriesExperiment::plotEnrichment(enrich=goa.fto, n_max = 10)
TimeSeriesExperiment::plotEnrichment(enrich=goa.pcdh, n_max = 10)

#---------------------------------------------------------------------------------------------------------------------------------------



