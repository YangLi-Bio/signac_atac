############################################################################
#                                                                          #
#            Filter differentiallt accessible regions in eGRNs             #
#                                                                          #
############################################################################


# 
# Run : Rscript Signac_regulons.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac/example.RDS hg38 ./
# 
# start_time <- Sys.time() # get the start time
message ("----> Loading libraries ...\n")
library(pbapply)
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/library.R") # libraries
source("/fs/ess/PCON0022/liyang/STREAM/param_tuning/codes/input.R") # Load data
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/globalVar.R") # global variables
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/functions.R") # define functions
# source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/filter_dorothea.R")
# filter eGRNs using the curated regulons in DorothEA

getGlobalVar() # build global variables
set.seed(1234) # set the seed for random computing
top.dars <- 100
top.egrns <- 1000
message ("----> Finished loading libraries.\n")


############################################################################
#                                                                          #
#                      Function to filter each eGRN                        #
#                                                                          #
############################################################################


filter_one_param <- function(data.name, param.name, org = "hg38", dbs = "KEGG_2019_Human") {
  
  # data.name : name of a dataset, e.g., 10x_pbmc_granulocyte_sorted_3k
  # param.name : name of a parameter setting, e.g., reso_0.8_fc_0.25_padj_0.01_var_0.25.qsave
  # org : the organism, e.g., hg38
  # dbs : the database for enrichment analysis, e.g., KEGG_2019_Human
  
  
  pbmc <- qs::qread(paste0(gsub(".qsave", "", param.name), "/", 
                           data.name, "_dir/", "pbmc_clusters.qsave"))
  out.dir <- paste0(gsub(".qsave", "", param.name), "/", 
                    data.name, "_dir/")
  
  
  # egrn.file <- paste0(data.name, "/", param.name) # paste the full name
  da_peaks <- qs::qread(paste0(gsub(".qsave", "", param.name), "/", 
                               data.name, "_dir/", "DARs.qsave"))
  
  
  # Annotate peaks using TF binding sites in JASPAR
  message ("----> Annotating peaks using JASPAR TF binding sites ...\n")
  hg.mus.jaspar.list <- qs::qread('/fs/ess/PCON0022/liyang/STREAM/databases/hg_mus_JASPAR_TF_binding_sites.qsave')
  # load the JASPAR TF binding sites
  
  double.hash <- add_TF(peaks = da_peaks$feature, hg.mus.jaspar.list = hg.mus.jaspar.list, 
                        org = org) # add TF binding sites
  peak.TFs <- double.hash$peak.hash # TFs binding each peak
  TF.peaks <- double.hash$TF.hash # peaks bound by each TF
  rm(hg.mus.jaspar.list)
  qs::qsave(double.hash, paste0(out.dir, "double_hash.qsave"))
  
  
  # Identify the DARs overlapped with gene promoters or distal cis-regulatory elements
  message ("----> Began identifying the DARs overlapped with gene promoters or distal cis-regulatory elements ...\n")
  # conns.use <- conns.clean[unlist(pbmclapply(1:nrow(conns.clean), function(i) {
  #   if(conns.clean$Peak1 %in% input.df$site_name | conns.clean$Peak2 %in% input.df$site_name) {
  #     return(T)
  #   } else {
  #     return(F)
  #   }
  # }, mc.cores = min(detectCores(), nrow(conns.clean))))] # only retain the peaks potentially linked to genes
  
  # distal.peaks <- setdiff(names(peak.TFs), input.df$site_name) # select the peaks not directly linked to genes
  # left.conns <- conns.clean[, c("Peak1", "Peak2")]
  # colnames(left.conns) <- c("distal", "site_name") # rename the columns
  # distal.conns <- left_join(x = left.conns, y = input.df, by = "site_name") %>% 
  #   dplyr::filter(!is.na(gene_name)) # get the peaks distally linked to genes
  # rm(left.conns)
  # distal.genes <- rbind(distal.conns[, c(1, 3)] %>% setNames(c("site_name", "gene_name")), input.df) %>% 
  #   split(f = .$site_name) %>% lapply(., "[[", ("gene_name")) # the relations between distal peaks and genes
  # qs::qsave(distal.genes, paste0(out.dir, "distal_genes.qsave"))
  # message ("----> Finished identified ", length(distal.genes), " pairs of peak-gene cis-regulatory relations.\n")
  distal.genes <- qs::qread(paste0(gsub(".qsave", "", param.name), "/", 
                            data.name, "_dir/", "distal_genes.qsave"))


  # Identify cell-type-active regulons
  message ("----> Began identifying cell-type-active regulons ...\n")
  regulons <- do.call(c, pbmclapply(unique(da_peaks$group), function(i) {
    dars.df <- da_peaks[da_peaks$group == i, , drop = F] # DARs in the cell cluster with at least one binding TFs
    dars <- dars.df %>% pull(feature)
    
  
    # Traverse all the TFs
    lapply(names(TF.peaks), function(j) {
      TF.dars <- intersect(dars, TF.peaks[[j]]) # DARs bound by this TF
      TF.dars.df <- dars.df[which(dars.df$feature %in% TF.dars),, drop = F]
      TF.dars.df <- TF.dars.df[order(TF.dars.df$logFC, decreasing = T),, drop = F]
      if (nrow(TF.dars.df) > top.dars) {
        TF.dars.df <- TF.dars.df[1:top.dars,, drop = F]
        TF.dars <- unique(TF.dars.df$feature)
      }
      return(list(TF = j, genes = Reduce(union, distal.genes[TF.dars]), 
                  peaks = TF.dars, cells = names(pbmc$orig.ident)[which(pbmc$seurat_clusters == i)]))
      # genes regulated by this TF
    })
  }, mc.cores = min(detectCores(), length(unique(da_peaks$group)))))
  
  size.ll <- sapply(regulons, function(x) {
    length(x$genes)
  })
  # names(size.ll) <- seq_along(regulons)
  order.ll <- order(size.ll, decreasing = T)
  regulons <- regulons[order.ll]
  regulons <- lapply(regulons, function(x) {
    if (length(x$genes) < 1 | length(x$peaks) < 1 | length(x$cells) < 1) {
      return(NULL)
    }
    x
  })
  regulons <- regulons[!sapply(regulons, is.null)]
  qs::qsave(regulons, paste0(out.dir, "regulons.qsave"))
  message ("----> Finished writing ", length(regulons), " eGRNs to ", out.dir, ".\n")
  
  
  # end_time <- Sys.time()
  # run_time <- end_time - start_time
  # message ("----> Running time: ", run_time, " min.\n")
  
  
  # Pathway enrichment analysis
  genes.ll <- lapply(regulons, "[[", "genes")
  pathway.ll <- Reduce(rbind, pblapply(seq_along(genes.ll), function(j) {
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
  }))
  # qs::qsave(pathway.ll, gsub(".qsave", ".extended_kegg.qsave", egrn.file))
  qs::qsave(pathway.ll, paste0(out.dir, "kegg_pathways.qsave"))
  
  
  # p.cutoff <- pval / n.pathways / length(genes.ll)
  # p.cutoff
  # kegg.df <- pathway.ll[pathway.ll$P.value < p.cutoff,]
  # precision <- length(unique(kegg.df$Regulon)) / length(genes.ll)
  # recall <- length(unique(kegg.df$Term)) / n.pathways
  # F1 <- 2 * precision * recall / (precision + recall)
  # 
  # length(genes.ll)
  # length(unique(kegg.df$Regulon))
  # length(unique(kegg.df$Term))
  # 
  # precision
  # recall
  # F1
  
}


############################################################################
#                                                                          #
#                                Main program                              #
#                                                                          #
############################################################################


# Get parameters
args <- commandArgs(T)
data.name <- args[1] # the name of a dataset, e.g., 10x_pbmc_granulocyte_sorted_3k
org <- args[2] # name of an organism, e.g., hg38


# Annotation databases
org.pathway <- orgPathway.hash[[org]] # organism ID in KEGG
org.gs <- orgGS.hash[[org]] # genome sequences
org.anno <- orgAnno.hash[[org]] # the annotation
org.db <- orgDB.hash[[org]] # database


# Convert characters into objects
specLib(org, org.anno, org.gs, org.db) # load specific libraries
org.gs <- get(org.gs)
org.anno <- get(org.anno)
org.db <- get(org.db)


# dbs <- "KEGG_2019_Mouse"
# dbs <- "KEGG_2019_Human"
if (grepl("^mm", org)) {
  dbs <- "KEGG_2019_Mouse"
  n.pathways <- 303
} else {
  dbs <- "KEGG_2019_Human"
  n.pathways <- 308
}


# param.ll <- list.files(path = data.name, pattern = ".kegg.csv") %>% 
#   gsub(pattern = ".kegg.csv", ".qsave", .) # get all files
param.ll <- setdiff(list.dirs(path = ".", recursive = F), 
                    c("./backup", "./Signac_atac")) %>% 
  gsub(pattern = "^\\./", "", .) %>% paste0(., ".qsave")
lapply(param.ll, function(x) {
  message ("----> Dataset: ", x, " ...\n")
  filter_one_param(data.name = data.name, param.name = x, org = org, dbs = dbs) # process one parameter setting
})
