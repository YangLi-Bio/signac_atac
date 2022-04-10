############################################################################
#                                                                          #
#          Filter eGRNs using the curated regulons in DoRothEA             #
#                                                                          #
############################################################################


# Load libraries
library(dorothea)
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/library.R")


# Filter eGRNs based on DoRothEA
# In this way, Signac-ATAC select CREs and genes based on both curated TF binding sites and regulons
filter_egrns_dorothea <- function(egrn.ll, org = "hg38") {
  
  # egrn.ll : the list of eGRNs, 
  #           e.g., /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/10x_pbmc_granulocyte_sorted_3k/reso_0.8_fc_0_padj_0.01_var_0.25.qsave
  # rds : the Seurat object file, e.g., 10x_pbmc_granulocyte_sorted_3k.RDS
  # org : the organism, e.g., hg38
  
  
  # Set the DoRothEA
  ifelse(grepl("^mm", org), db <- "dorothea_mm", db <- "dorothea_hs")
  regulon.ll <- get(db) %>% 
    filter(confidence %in% c("A", "B")) %>% split(., f = .$tf) %>% sapply(., "[[", "target")
  filtered.egrns <- pbmclapply(egrn.ll, function(x) {
    if (!x$TF %in% names(regulon.ll)) {
      return(NULL)
    }
    genes <- intersect(x$genes, regulon.ll[[x$TF]]) # filter genes using DoRothEA regulons
    x$genes <- genes # update the genes
    x
  }, mc.cores = detectCores())
  filtered.egrns <- filtered.egrns[!sapply(filtered.egrns, function(x) {
    is.null(x) | length(x$genes) < 1
  })] # remove NULL
  egrn.size <- sapply(filtered.egrns, function(x) {
    return(length(x$genes))
  }) # count the number of genes in each eGRN
  filtered.egrns <- filtered.egrns[order(egrn.size, decreasing = T)]
}


# Enrichment analysis
# dbs <- "KEGG_2019_Human"
genes.ll <- lapply(filtered.egrns, "[[", "genes")
pathway.ll <- Reduce(rbind, pbmclapply(seq_along(genes.ll), function(j) {
  message ("------------> Performing enrichment analysis for regulon: ", j, " ...\n")
  genes <- genes.ll[[j]]
  
  if (length(genes) < 1) {
    return(NULL)
  }
  
  pathways <- tryCatch(enrichr(genes, dbs) %>% `[[` (dbs), 
                       error = function(e) {
                         0 })
  
  if (is.null(pathways) | length(pathways) < 1 | is.numeric(pathways) | 
      nrow(pathways) < 1 | ncol(pathways) != 9) {
    return(NULL)
  }
  
  cbind(Regulon = rep(j, nrow(pathways)), pathways)
}, mc.cores = detectCores()))
