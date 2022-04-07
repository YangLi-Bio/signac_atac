################################################################################
#                                                                              #
#            Generate KEGG pathway enrichment analysis results                 #
#                                                                              #
################################################################################


# Load libraries
library(pbmcapply)
library(pbapply)
library(qs)
library(enrichR)
library(dplyr)


# Global variables
args <- commandArgs(T)
dd <- args[1]
org <- args[2]

root.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac"
datalist.file <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt"
message ("----> Root directory: ", root.dir, ".\n")


message ("--------> Processing the dataset ", dd, " ...\n")
file.ll <- Sys.glob(file.path(dd, "*.qsave")) # regulons for different parameter settings

if (length(file.ll) < 1) {
  next
}

dbs <- ifelse(grepl("^mm", org), "KEGG_2019_Mouse", "KEGG_2019_Human")

pblapply(file.ll, function(ff) {
  regulon.ll <- qs::qread(ff) # regulon list
  flag.ll <- sapply(lapply(regulon.ll, "[[", "genes"), length) > 0
  regulon.ll <- regulon.ll[flag.ll] # filter regulons
  genes.ll <- lapply(regulon.ll, "[[", "genes")
  # list of regulon genes for one dataset
  
  if (length(genes.ll) < 1) {
    return(NULL)
  }
  
  pathway.ll <- Reduce(rbind, pbmclapply(seq_along(genes.ll), function(j) {
    message ("------------> Performing enrichment analysis for regulon: ", j, " ...\n")
    genes <- genes.ll[[j]]
    
    if (length(genes) < 1) {
      return(NULL)
    }
    
    pathways <- tryCatch(enrichr(genes, dbs) %>% `[[` (dbs), 
                         error = function(e) {
                           0 })
    
    if (is.null(pathways) | is.numeric(pathways) | nrow(pathways) == 0) {
      return(NULL)
    }
    
    cbind(Regulon = rep(j, nrow(pathways)), pathways)
  }, mc.cores = min(length(file.ll), detectCores())))
  
  write.csv(pathway.ll[, c(-5, -6)], gsub(".qsave", ".kegg.csv", ff))
})