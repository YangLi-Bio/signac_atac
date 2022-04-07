############################################################################
#                                                                          #
#                           Load the R libraries                           #
#                                                                          #
############################################################################

library(cicero) # co-accessibility analysis
message ("cicero\n")

library(Seurat) # single-cell RNA-Seq and multimodal analysis tools
library(Signac) # single-cell ATAC-Seq or multimodal analysis tools
library(GenomeInfoDb) # the genome information annotation
library(ggplot2) # for plotting
library(patchwork) # parallel computing
library(hash) # hash table
library(dplyr) # pipe functions
message ("Finished basic libraries\n")

library(pbmcapply) # parallel computing
library(qs) # fast loading and saving
library(Repitools) # epigenetic tools
library(GenomicRanges) # genomic ranges
library(data.table) # convert nested list into data frame
# library(SeuratWrappers) # a collection of community-provided methods and extensions for Seurat


# add specific libraries
specLib <- function(org = 'hg38', org.anno = 'EnsDb.Hsapiens.v86', 
                    org.gs = 'BSgenome.Hsapiens.UCSC.hg38', 
                    org.db = 'org.Hs.eg.db') {
  library(org.anno, character.only = T)
  library(org.gs, character.only = T)
  library(org.db, character.only = T)
}

message ("Finished all.\n")