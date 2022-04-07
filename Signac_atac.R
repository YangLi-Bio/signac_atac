############################################################################
#                                                                          #
#                           Load the R codes                               #
#                                                                          #
############################################################################


# 
# Run : Rscript Signac_regulons.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac/example.RDS hg38 ./
# 
start_time <- Sys.time() # get the start time
message ("----> Loading libraries ...\n")
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/library.R") # libraries
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/globalVar.R") # global variables
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/functions.R") # define functions
source("/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac_atac/filter_dorothea.R")
# filter eGRNs using the curated regulons in DorothEA

getGlobalVar() # build global variables
set.seed(1234) # set the seed for random computing
message ("----> Finished loading libraries.\n")


############################################################################
#                                                                          #
#                           Load the dataset                               #
#                                                                          #
############################################################################

# 
# Run : Rscript Signac_atac.R /fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac/example.RDS hg38 ./
# 

# The format of result directory: reso_0.8_fc_0_padj_0.01_var_0
args <- commandArgs(T) # get parameters from Shell command
rds.path <- args[1] # the path to RDS file, e.g, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/10x_human_brain_3k.RDS"
org <- args[2] # the organism, e.g., "hg38"
out.dir <- args[3] # the output directory

# Parameters to tune
clust.resolution <- as.numeric(args[4]) # the resolution for clustering
log.fc.cutoff <- as.numeric(args[5]) # the cutoff of log fold change when identifying DARs
padj.cutoff <- as.numeric(args[6]) # the cutoff of adjusted p-values when predicting DARs
covar.cutoff <- as.numeric(args[7]) # the cutoff of covariance when running Cicero


# rds.path <- "/fs/ess/PCON0022/liyang/STREAM/benchmarking/Signac/example.RDS"
# org <- "hg38"
# out.dir <- "./"


out.dir.array <- strsplit(x = rds.path, split = "[/|\\.]") %>% `[[` (1) # split the input file
# message ("Length: ", length(out.dir.array), "\n")
out.dir <- paste0(out.dir, out.dir.array[length(out.dir.array) - 1], "_dir/")
dir.create(path = out.dir)
# create the output directory


org.pathway <- orgPathway.hash[[org]] # organism ID in KEGG
org.gs <- orgGS.hash[[org]] # genome sequences
org.anno <- orgAnno.hash[[org]] # the annotation
org.db <- orgDB.hash[[org]] # database


# convert characters into objects
specLib(org, org.anno, org.gs, org.db) # load specific libraries
org.gs <- get(org.gs)
org.anno <- get(org.anno)
org.db <- get(org.db)


pbmc <- readRDS(file = rds.path) # read the RDS file
message ("----> Finished loading the dataset composed of ", ncol(pbmc), " cells.\n")
DefaultAssay(pbmc) <- "ATAC"

# Remove the RNA assay
if(!is.null(pbmc@assays$RNA)) {
  pbmc[["RNA"]] <- NULL # only retain the "ATAC" assay
}

# pbmc <- subset(x = pbmc, features = rownames(pbmc)[1:5000])


############################################################################
#                                                                          #
#                   Dimension reduction and clustering                     #
#                                                                          #
############################################################################


# dimension reduction
message ("----> Performing dimension reduction ...\n")
pbmc <- RunTFIDF(pbmc) # normalizarion
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # retain the top features according to quantiles
pbmc <- RunSVD(pbmc) # singular value decomposition

# correlation with sequencing depth
png(filename = paste0(out.dir, "depth_corr.png"))
DepthCor(pbmc)
dev.off()

# clustering
message ("----> Performing cell clustering ...\n")
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = clust.resolution)
png(filename = paste0(out.dir, "clusters.png"))
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()
qs::qsave(pbmc, paste0(out.dir, "pbmc_clusters.qsave"))
message ("----> Finished identifying ", length(unique(pbmc$seurat_clusters)), " cell types.\n")


############################################################################
#                                                                          #
#                      Co-accessibility analysis                           #
#                                                                          #
############################################################################


message ("----> Began running Cicero to build peak co-accessibility linkages.\n")
m.atac <- GetAssayData(object = pbmc, slot = "data", assay = "ATAC") # get the accessibility matrix

# Make CDS data
summ <- summary(m.atac) # convert the matrix into a sparse matrix
cicero.data <- data.frame(Origin = rownames(m.atac)[summ$i],
                          Destination = colnames(m.atac)[summ$j],
                          Weight      = summ$x) # transform the sparse matrix into a data frame
input.cds <- make_atac_cds(cicero.data, binarize = F) %>% detect_genes
# build the CDS data based on the data frame


# Run Cicero
#
# Since Cicero failed in generating co-accessibility links on some small cell clusters, we ran Cicero
# Cicero on the whole ATAC matrix
#
input.cds <- input.cds[Matrix::rowSums(exprs(input.cds)) != 0, ] %>% estimate_size_factors %>% 
  preprocess_cds(method = "LSI", verbose = F) %>% reduce_dimension(reduction_method = 'UMAP', 
                                                                   preprocess_method = "LSI")
qs::qsave(input.cds, paste0(out.dir, "input_cds.qsave"))
# input.cds <- estimate_size_factors(input.cds) # estimate size factor
# input.cds <- preprocess_cds(input.cds, method = "LSI", verbose = F) # data transformation
# input.cds <- reduce_dimension(input.cds, reduction_method = 'UMAP', 
#                               preprocess_method = "LSI") # dimensional reduction

umap.coords <- reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
cicero.cds <- make_cicero_cds(input.cds, reduced_coordinates = umap.coords)
genome.info <- data.frame(org.gs@seqinfo@seqnames, 
                          org.gs@seqinfo@seqlengths) # genome sequence lengths
colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
# genome.info$seqlengths <- 10000000
# genome.info <- subset(genome.info, seqnames == "chr1") # for debugging
# conns <- run_cicero(cicero.cds, genome.info, sample_num = 2, 
#                     window = 500000) # build peak-peak linkages using cicero
conns <- run_cicero(cicero.cds, genome.info,
                    window = 500000) # build peak-peak linkages using cicero
# ccans <- generate_ccans(conns) # generate CCANs
# Some redundant co-accessibility linkages exist

conns$Peak2 <- as.character(conns$Peak2) # unify the format of columns
conns$coaccess <- as.numeric(conns$coaccess) # convert characters into numeric values
qs::qsave(conns, paste0(out.dir, "conns.qsave")) # save the co-accessibility linkages
message ("----> Identified ", nrow(conns), " co-accessibility linkages.\n")


# Rearrange the peaks in each row according to alphabetical order
message ("----> Began removing duplicated co-accessibility linkages ...\n")
conns.clean <- conns[conns$coaccess > covar.cutoff,]
# conns.clean <- rbindlist(pbmclapply(1:nrow(conns), function(i) {
#   quiet(ifelse(conns$Peak1[i] >= conns$Peak2[i], 
#          rr <- list(Peak1 = conns$Peak1[i], 
#                  Peak2 = conns$Peak2[i], 
#                  coaccess = conns$coaccess[i]), 
#          rr <- list(Peak1 = conns$Peak2[i], 
#                  Peak2 = conns$Peak1[i], 
#                  coaccess = conns$coaccess[i]))) # reorder the columns
#   
#   return(rr)
# }, mc.cores = min(detectCores(), nrow(conns)))) %>% distinct
# qs::qsave(conns.clean, paste0(out.dir, "conns_clean.qsave")) # save the co-accessibility linkages

conns.clean <- conns
rm(m.atac)
rm(conns)
message ("----> Retained ", nrow(conns.clean), " co-accessibility linkages.\n")


# An example
# conns <- run_cicero(cicero.cds, subset(genome.info, seqnames == "chr1"),
#                     window = 500000, sample_num = 2)

# Links(pbmc) <- ConnectionsToLinks(conns = conns) # add links to Seurat Object


# Run Cicero in each subset of cells
# # Build co-accessibility connections
# pbmclapply(seq_along(levels(pbmc$seurat_clusters)), function(i) {
#   id <- as.numeric(levels(pbmc$seurat_clusters)[i]) # cluster ID
#   m.atac <- big.m[, colnames(pbmc)[pbmc$seurat_clusters == id]] # get the sub-matrices
#   
#   
#   # Make CDS data
#   summ <- summary(m.atac) # convert the matrix into a sparse matrix
#   cicero.data <- data.frame(Origin = rownames(m.atac)[summ$i],
#                             Destination = colnames(m.atac)[summ$j],
#                             Weight      = summ$x) # transform the sparse matrix into a data frame
#   input.cds <- make_atac_cds(cicero.data, binarize = F) %>% detect_genes
#   # build the CDS data based on the data frame
# 
#   
#   # Run Cicero
#   # input.cds <- detect_genes(input.cds) # set the global expression detection threshold
#   input.cds <- input.cds[Matrix::rowSums(exprs(input.cds)) != 0, ] %>% estimate_size_factors %>% 
#     preprocess_cds(method = "LSI", verbose = F) %>% reduce_dimension(reduction_method = 'UMAP', 
#                                                                     preprocess_method = "LSI")
# 
#   # input.cds <- estimate_size_factors(input.cds) # estimate size factor
#   # input.cds <- preprocess_cds(input.cds, method = "LSI", verbose = F) # data transformation
#   # input.cds <- reduce_dimension(input.cds, reduction_method = 'UMAP', 
#   #                               preprocess_method = "LSI") # dimensional reduction
#   
#   umap.coords <- reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
#   cicero.cds <- make_cicero_cds(input.cds, reduced_coordinates = umap.coords)
#   genome.info <- data.frame(org.gs@seqinfo@seqnames, 
#                             org.gs@seqinfo@seqlengths) # genome sequence lengths
#   colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
#   conns <- run_cicero(cicero.cds, genome.info,
#                       window = 500000) # build peak-peak linkages using cicero
#  
#   # conns <- run_cicero(cicero.cds, subset(genome.info, seqnames == "chr1"),
#   #                     window = 500000, sample_num = 2)
#    
# }, mc.cores = min(detectCores(), length(levels(pbmc$seurat_clusters))))


############################################################################
#                                                                          #
#                      Build peak-gene linkages                            #
#                                                                          #
############################################################################


message ("----> Began building peak-gene linkages ...\n")
exon.anno <- Annotation(pbmc) # the exon annotations
gene_anno <- gr2dt(exon.anno) # convert GRanges into data frame
# names(exon.anno) <- NULL # avoid duplicated names
# gene_anno <- annoGR2DF(exon.anno) # get the gene information


#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 

# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$tx_id),]

# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 

# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$tx_id),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[, c("seqnames", "start", "end", "gene_name")]
input.cds <- annotate_cds_by_site(input.cds, gene_annotation_sub)
input.df <- fData(input.cds) %>% as.data.frame() %>% dplyr::filter(!is.na(gene_name)) %>% 
  dplyr::select("site_name", "gene_name")
# get DFrame to save the connections between peaks and genes

qs::qsave(input.df, paste0(out.dir, "input_df.qsave"))
message ("----> Finished building ", nrow(input.df), " peak-gene linkages.\n")


############################################################################
#                                                                          #
#                   Differential accessibility analysis                    #
#                                                                          #
############################################################################


message ("----> Begain performing differential accessibility analysis ...\n")
da_peaks <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'seurat_clusters', 
                                          assay = 'data', seurat_assay = 'ATAC') %>% 
  # dplyr::filter(padj <= 0.05, logFC >= 0.25, pct_in >= 0.25, pct_out >= 0.25)
  dplyr::filter(padj <= padj.cutoff, logFC > log.fc.cutoff)
# find differentially accessible regions (DARs)

# TBD : maybe use a stringent standards to filter DARs!
qs::qsave(da_peaks, paste0(out.dir, "DARs.qsave"))
message ("----> Discovered ", nrow(da_peaks), " differentially accessible regions (DARs).\n")


# Annotate peaks using TF binding sites in JASPAR
message ("----> Annotating peaks using JASPAR TF binding sites ...\n")
hg.mus.jaspar.list <- qs::qread('/fs/ess/PCON0022/liyang/STREAM/databases/hg_mus_JASPAR_TF_binding_sites.qsave')
# load the JASPAR TF binding sites

double.hash <- add_TF(peaks = da_peaks$feature, hg.mus.jaspar.list = hg.mus.jaspar.list, 
                      org = org) # add TF binding sites
peak.TFs <- double.hash$peak.hash # TFs binding each peak
TF.peaks <- double.hash$TF.hash # peaks bound by each TF
rm(hg.mus.jaspar.list)


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
left.conns <- conns.clean[, c("Peak1", "Peak2")]
colnames(left.conns) <- c("distal", "site_name") # rename the columns
distal.conns <- left_join(x = left.conns, y = input.df, by = "site_name") %>% 
  dplyr::filter(!is.na(gene_name)) # get the peaks distally linked to genes
rm(left.conns)
distal.genes <- rbind(distal.conns[, c(1, 3)] %>% setNames(c("site_name", "gene_name")), input.df) %>% 
  split(f = .$site_name) %>% lapply(., "[[", ("gene_name")) # the relations between distal peaks and genes
qs::qsave(distal.genes, paste0(out.dir, "distal_genes.qsave"))
message ("----> Finished identified ", length(distal.genes), " pairs of peak-gene cis-regulatory relations.\n")


# Identify cell-type-active regulons
message ("----> Began identifying cell-type-active regulons ...\n")
regulons <- do.call(c, pbmclapply(unique(da_peaks$group), function(i) {
  dars <- da_peaks[da_peaks$group == i, ] %>% pull(feature) # DARs in the cell cluster with at least one binding TFs
  
  # Traverse all the TFs
  lapply(names(TF.peaks), function(j) {
    TF.dars <- intersect(dars, TF.peaks[[j]]) # DARs bound by this TF
    return(list(TF = j, genes = Reduce(union, distal.genes[TF.dars]), 
                peaks = TF.dars, cells = names(pbmc$orig.ident)[which(pbmc$seurat_clusters == i)]))
    # genes regulated by this TF
  })
}, mc.cores = min(detectCores(), length(unique(da_peaks$group)))))
qs::qsave(regulons, paste0(out.dir, "regulons.qsave"))
message ("----> Finished identifying ", length(regulons), " cell-type-specific regulons.\n")


end_time <- Sys.time()
run_time <- end_time - start_time
message ("----> Running time: ", run_time, " min.\n")