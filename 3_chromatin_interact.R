################################################################################
#                                                                              #
#            Generate enhancer-target interaction-based results                #
#                                                                              #
################################################################################


# Load libraries
library(pbmcapply)
library(qs)
library(dplyr)
library(data.table)
library(tidyverse)


# Global variables
root.dir <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac"
#datalist.file <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/dataset_list.txt" # list of datasets
interact.file <- "/fs/ess/PCON0022/liyang/STREAM/comparison/EnhancerAtlas2/enhancer_target_interact.txt" # db file
db.path <- "/fs/ess/PCON0022/liyang/STREAM/comparison/EnhancerAtlas2/curated_interactions"
# path of the curated enhancer-target interactions

message ("----> Root directory: ", root.dir, ".\n")


# Load the result directories
data.org.db.ll <- strsplit(readLines(interact.file), split = " ") # the file to save dataset, organism, and db files
result.dir <- paste(root.dir, lapply(data.org.db.ll, "[[", 1), sep = "/") # result directories
#org.ll <- unlist(lapply(data.org.db.ll, "[[", 2)) # organism
interact.db.ll <- unlist(lapply(data.org.db.ll, "[[", 3)) # db files
message ("----> There are results on ", length(result.dir), " datasets.\n")


# Calculate the ratio of enhancer-target interactions incorporated in CREs underlying regulons
message ("----> Began calculating the ratio of enhancer-target interactions ...\n")
pblapply(seq_along(result.dir), function(i) {
  db <- interact.db.ll[[i]]
  
  if (db == "-") {
      return(NULL)
  }
  
  dd <- result.dir[[i]]
  message ("--------> Processing the dataset ", dd, " ...\n")
  file.ll <- Sys.glob(file.path(dd, "*.qsave")) # regulons for different parameter settings
  
  if (length(file.ll) < 1) {
    next
  }
  
  
  # Load the curated enhancer-target interactions
  message ("--------> Loading the enhancer-target databases ...\n")
  db.ll <- strsplit(db, ",")[[1]] # in case of multiple files
  db.interact <- rbindlist(do.call("c", lapply(db.ll, function(x) {
    ln.ll <- read.table(paste0(db.path, "/", x)) %>% pull(V1) %>% strsplit(., split = "_")
    return(lapply(ln.ll, function(y) {
      peak <- gsub(pattern = "\\:", replacement = "-", y[[1]])
      gene <- strsplit(y[[2]], split = "\\$") %>% "[[" (1) %>% "[[" (2)
      return(list(peak = peak, gene = gene))
        }))
  }))) %>% distinct
  
  
  lapply(file.ll, function(ff) {
    regulon.ll <- qs::qread(ff) # regulon list
    flag.ll <- sapply(lapply(regulon.ll, "[[", "genes"), length) > 0
    regulon.ll <- regulon.ll[flag.ll] # filter regulons
    genes.ll <- lapply(regulon.ll, "[[", "genes")
    # list of regulon genes for one dataset
    
    peaks.ll <- lapply(regulon.ll, "[[", "peaks")
    # list of regulon peaks for one dataset
    
    if (length(genes.ll) < 1) {
      return(NULL)
    }
    
    pathway.ll <- Reduce(rbind, lapply(seq_along(genes.ll), function(j) {
      genes <- genes.ll[[j]]
      
      if (length(genes) < 1) {
        return(NULL)
      }
      
      pathways <- enrichr(genes, dbs) %>% `[[` (dbs)
      cbind(Regulon = rep(j, nrow(pathways)), pathways)
    }))
    
    write.csv(pathway.ll[, c(-5, -6)], gsub(".qsave", ".kegg.csv", ff))
  })
})
message ("----> Finished performing pathway enrichment analyses on ", length(result.dir), " datasets.\n")