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
root.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac"
datalist.file <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt"
message ("----> Root directory: ", root.dir, ".\n")


# Load the result directories
# result.dir <- list.dirs(root.dir) %>% 
#   setdiff(paste(root.dir, "scripts", sep = "/")) %>% 
#   setdiff(root.dir)

datalist.ar <- strsplit(readLines(datalist.file), split = " ")
result.dir <- paste(root.dir, lapply(datalist.ar, "[[", 1), sep = "/") # result directories
org.ll <- unlist(lapply(datalist.ar, "[[", 2)) # organism
message ("----> There are results on ", length(result.dir), " datasets.\n")


# KEGG enrichment analyses
message ("----> Began performing enrichment analyses on different parameter settings for multiple datasets ...\n")
pbmclapply(1:length(result.dir), function(i) {
  dd <- result.dir[[i]]
  message ("--------> Processing the dataset ", dd, " ...\n")
  file.ll <- Sys.glob(file.path(dd, "*.qsave")) # regulons for different parameter settings
  
  if (length(file.ll) < 1) {
    next
  }
  
  dbs <- ifelse(grepl("^mm", org.ll[[i]]), "KEGG_2019_Mouse", "KEGG_2019_Human")
  
  lapply(file.ll, function(ff) {
    regulon.ll <- qs::qread(ff) # regulon list
    flag.ll <- sapply(lapply(regulon.ll, "[[", "genes"), length) > 0
    regulon.ll <- regulon.ll[flag.ll] # filter regulons
    genes.ll <- lapply(regulon.ll, "[[", "genes")
    # list of regulon genes for one dataset
    
    if (length(genes.ll) < 1) {
      return(NULL)
    }
    
    pathway.ll <- Reduce(rbind, lapply(seq_along(genes.ll), function(j) {
      message ("------------> Performing enrichment analysis for regulon: ", j, " ...\n")
      genes <- genes.ll[[j]]
      
      if (length(genes) < 1) {
        return(NULL)
      }
      
      pathways <- enrichr(genes, dbs) %>% `[[` (dbs)
      cbind(Regulon = rep(j, nrow(pathways)), pathways)
    }))
    
    write.csv(pathway.ll[, c(-5, -6)], gsub(".qsave", ".kegg.csv", ff))
  })
}, mc.cores = min(length(result.dir), detectCores()))
message ("----> Finished performing pathway enrichment analyses on ", length(result.dir), " datasets.\n")