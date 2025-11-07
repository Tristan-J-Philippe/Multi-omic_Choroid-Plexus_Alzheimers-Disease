###
### 081 Spatial_cellchat.R
### This script was written by Dr. Denis R Avey
# Purpose: Calculate spatial proximity and distance-based features between cell types.
# Dependencies:
source("/080 Spatial_cellchat_util.R", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder



#Script of Spatial cell chat for Niche1 AD, run for all niches and groups.
# Load the human ligand-receptor database
data(CellChatDB.human)

# Define the desired order
new_levels <- c("Endothelial", "Epithelial", "Fibroblast", "Mural", "Immune")

# load Seurat object and ensure cell type levels are in the desired order
ChP_AD_niche1 <- readRDS("ChP_AD_niche1.rds")
ChP_AD_niche1$newMajorcelltype <- droplevels(ChP_AD_niche1$newMajorcelltype)
ChP_AD_niche1$newMajorcelltype <- factor(ChP_AD_niche1$newMajorcelltype, levels = new_levels)

# convert ChP metadata to matrix
meta_data <- ChP_AD_niche1@meta.data  

# extract x-y coordinates, subset by AD/AD
xy1 <- ChP_AD_niche1@images[["Slide.1"]]@boundaries[["centroids"]]@coords
xy2 <- ChP_AD_niche1@images[["Slide.2"]]@boundaries[["centroids"]]@coords
xy <- rbind(xy1,xy2)
rownames(xy) <- meta_data$unique_cellID

# calculate spot.size for AD and AD; use same value (3.48) for each
conversion.factor = 0.18
d = computeCellDistance(xy)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
# minimum spot.size set to 3.48 after empirical investigation
spatial.factors = data.frame(ratio = conversion.factor, tol = 3.48/2)

# Function to create CellChat objects; SCT version
create_cellchat <- function(seurat_obj, metadata, coordinates, spatial_factors) {
  DefaultAssay(seurat_obj) <- "SCT_full"
  cellchat_obj <- createCellChat(object = seurat_obj, meta = metadata, 
                                 group.by = "newMajorcelltype", 
                                 datatype = "spatial", coordinates = coordinates,
                                 spatial.factors = spatial_factors, assay = "SCT_full")
  return(cellchat_obj)
}

# Test on single AD/AD individual
cellchat_AD <- create_cellchat(ChP_AD_niche1, meta_data, xy, spatial.factors)
cellchat_AD@meta[["samples"]] <- droplevels(cellchat_AD@meta[["samples"]])

# Add harmony embeddings to CellChat object
cellchat_AD@meta$Harmony1 <- ChP_AD_niche1@reductions[["full_harmony_umap"]]@cell.embeddings[,1]
cellchat_AD@meta$Harmony2 <- ChP_AD_niche1@reductions[["full_harmony_umap"]]@cell.embeddings[,2]

# remove ChP object and metadata to save RAM
rm(ChP_AD_niche1)
rm(meta_data)
rm(xy)
gc()

# load and visualize all Cellchat genes by category
CellChatDB <- CellChatDB.human
# test first on a smaller subset of L-R genes (e.g. ECM-Receptor, Secreted Signaling, or Cell-Cell Contact)
# CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor", key = "annotation") # use ECM-Receptor, single category selected

cellchat_AD@DB <- CellChatDB
cellchat_AD <- subsetData(cellchat_AD) # subsetting for all gene sets

# you need presto library to run these functions
options(future.globals.maxSize = 20*1024^3)
future::plan("multisession", workers = 6)
cellchat_AD <- identifyOverExpressedGenes(cellchat_AD)
cellchat_AD <- identifyOverExpressedInteractions(cellchat_AD)
saveRDS(cellchat_AD, file = "cellchat_AD_SCT_niche1_initial.rds") # save intermediate

# this step takes a lot of computational resources
cellchat_AD <- computeCommunProb(cellchat_AD, trim = 0.05, type = "truncatedMean", contact.range = 20, contact.knn.k = NULL, scale.distance = 0.69)
# complete analysis
cellchat_AD <- filterCommunication(cellchat_AD, min.cells = 10) # This filter requires minimum of 10 cells per cell type
df.net_AD <- subsetCommunication(cellchat_AD)
head(df.net_AD)
cellchat_AD <- computeCommunProbPathway(cellchat_AD)
cellchat_AD <- aggregateNet(cellchat_AD)

# save cellchat objects
saveRDS(cellchat_AD, file = "cellchat_AD_SCT_niche1.rds")
#repeat for other niches and groups.

# read in CellChat objects (unnormalized)
NCI <- readRDS(file="ChP_NCI_cellchat_SCT_raw.rds")
AD <- readRDS(file="ChP_AD_cellchat_SCT_raw.rds")
NCI_niche1 <- readRDS(file="cellchat_NCI_SCT_niche1.rds")
AD_niche1 <- readRDS(file="cellchat_AD_SCT_niche1.rds")
NCI_niche2 <- readRDS(file="cellchat_NCI_SCT_niche2.rds")
AD_niche2 <- readRDS(file="cellchat_AD_SCT_niche2.rds")
NCI_niche3 <- readRDS(file="cellchat_NCI_SCT_niche3.rds")
AD_niche3 <- readRDS(file="cellchat_AD_SCT_niche3.rds")

# run function (example shown for niche 1)
result_NCI1_pot <- ComputePotentialSpatialFractions(
  NCI_niche1,
  n_workers = 4,
  log_file = "NCI_k3_n1_potential_progress.log",
  output_csv = "NCI_k3_n1_potentialProportions.csv"
)

result_AD1_pot <- ComputePotentialSpatialFractions(
  AD_niche1,
  n_workers = 4,
  log_file = "AD_k3_n1_potential_progress.log",
  output_csv = "AD_k3_n1_potentialProportions.csv"
)
# repeat for other cellchat objects



# color palette for cell types
major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B") 

# load CC objects
NCI <- readRDS(file="data_objects/ChP_NCI_cellchat_SCT_raw.rds")
AD <- readRDS(file="data_objects/ChP_AD_cellchat_SCT_raw.rds")
NCI_niche1 <- readRDS(file="data_objects/ChP_NCI_cellchat_SCT_niche1.rds")
AD_niche1 <- readRDS(file="data_objects/ChP_AD_cellchat_SCT_niche1.rds")
NCI_niche2 <- readRDS(file="data_objects/ChP_NCI_cellchat_SCT_niche2.rds")
AD_niche2 <- readRDS(file="data_objects/ChP_AD_cellchat_SCT_niche2.rds")
NCI_niche3 <- readRDS(file="data_objects/ChP_NCI_cellchat_SCT_niche3.rds")
AD_niche3 <- readRDS(file="data_objects/ChP_AD_cellchat_SCT_niche3.rds")

# potential interactions given paired cell distances, agnostic to interaction_name
AD_potential <- read.csv("data_csv/AD_potentialProportions.csv")
NCI_potential <- read.csv("data_csv/NCI_potentialProportions.csv")
niche1_AD_potential <- read.csv("data_csv/AD_k3_n1_potentialProportions.csv")
niche1_NCI_potential <- read.csv("data_csv/NCI_k3_n1_potentialProportions.csv")
niche2_AD_potential <- read.csv("data_csv/AD_k3_n2_potentialProportions.csv")
niche2_NCI_potential <- read.csv("data_csv/NCI_k3_n2_potentialProportions.csv")
niche3_AD_potential <- read.csv("data_csv/AD_k3_n3_potentialProportions.csv")
niche3_NCI_potential <- read.csv("data_csv/NCI_k3_n3_potentialProportions.csv")


# Normalize CellChat object with spatially aware potential interaction counts
cellchat_nci_norm <- normalize_cellchat_by_log_npairs(NCI, NCI_potential) 
cellchat_ad_norm <- normalize_cellchat_by_log_npairs(AD, AD_potential)
cellchat_nci_niche1_norm <- normalize_cellchat_by_log_npairs(NCI_niche1, niche1_NCI_potential) 
cellchat_ad_niche1_norm <- normalize_cellchat_by_log_npairs(AD_niche1, niche1_AD_potential)
cellchat_nci_niche2_norm <- normalize_cellchat_by_log_npairs(NCI_niche2, niche2_NCI_potential)
cellchat_ad_niche2_norm <- normalize_cellchat_by_log_npairs(AD_niche2, niche2_AD_potential)
cellchat_nci_niche3_norm <- normalize_cellchat_by_log_npairs(NCI_niche3, niche3_NCI_potential)
cellchat_ad_niche3_norm <- normalize_cellchat_by_log_npairs(AD_niche3, niche3_AD_potential)

# save normalized objects
saveRDS(cellchat_nci_norm, file = "cellchat_nci_norm.rds")
saveRDS(cellchat_ad_norm, file = "cellchat_ad_norm.rds")
saveRDS(cellchat_nci_niche1_norm, file = "cellchat_nci_k3_niche1_norm.rds")
saveRDS(cellchat_ad_niche1_norm, file = "cellchat_ad_k3_niche1_norm.rds")
saveRDS(cellchat_nci_niche2_norm, file = "cellchat_nci_k3_niche2_norm.rds")
saveRDS(cellchat_ad_niche2_norm, file = "cellchat_ad_k3_niche2_norm.rds")
saveRDS(cellchat_nci_niche3_norm, file = "cellchat_nci_k3_niche3_norm.rds")
saveRDS(cellchat_ad_niche3_norm, file = "cellchat_ad_k3_niche3_norm.rds")

# Compare signaling strength and number of interactions

# Define the updated CellChat object list
cellchat_list <- list(
  NCI = NCI,
  AD = AD,
  NCI_n1 = NCI_niche1,
  AD_n1 = AD_niche1,
  NCI_n2 = NCI_niche2,
  AD_n2 = AD_niche2,
  NCI_n3 = NCI_niche3,
  AD_n3 = AD_niche3
)

norm_cellchat_list <- list(
  NCI = cellchat_nci_norm,
  AD = cellchat_ad_norm,
  NCI_n1 = cellchat_nci_niche1_norm,
  AD_n1 = cellchat_ad_niche1_norm,
  NCI_n2 = cellchat_nci_niche2_norm,
  AD_n2 = cellchat_ad_niche2_norm,
  NCI_n3 = cellchat_nci_niche3_norm,
  AD_n3 = cellchat_ad_niche3_norm
)

# Define condition names
conditions <- names(norm_cellchat_list)

# Compute overall signaling strength (normalized)
signal_strengths <- sapply(norm_cellchat_list, function(cellchat) {
  if (!is.null(cellchat@net$prob)) {
    return(sum(cellchat@net$prob, na.rm = TRUE))
  } else {
    return(NA)
  }
})

# Compute total number of interactions (unnormalized)
num_interactions <- sapply(cellchat_list, function(cellchat) {
  if (!is.null(cellchat@net$count)) {
    return(sum(cellchat@net$count, na.rm = TRUE))
  } else {
    return(NA)
  }
})

# AD-differential (AD - NCI) barplot of signaling strength 

# Compute AD-differential signaling strength for each pair
fc_data <- lapply(names(condition_pairs), function(name) {
  ad_val <- signal_strengths[condition_pairs[[name]][1]]
  nci_val <- signal_strengths[condition_pairs[[name]][2]]
  diff <- ad_val - nci_val
  data.frame(Comparison = name, Diff = diff)
})
df_fc <- do.call(rbind, fc_data)

# Define custom colors for comparisons
fc_colors <- c(
  Full      = "#E41A1C",     
  Niche1    = "#33a02c",   
  Niche2    = "#1f78b4",
  Niche3    = "#ffb000"
)

ggplot(df_fc, aes(x = Comparison, y = Diff, fill = Comparison)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = fc_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "AD-differential Normalized Signaling Strength",
    x = NULL,
    y = "(AD - NCI)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")

# merge and save cellchat objects
cellchat_merged <- mergeCellChat(list(NCI, AD), add.names = c("NCI","AD"))
cellchat_merged1 <- mergeCellChat(list(NCI_niche1, AD_niche1), add.names = c("NCI_1","AD_1"))
cellchat_merged2 <- mergeCellChat(list(NCI_niche2, AD_niche2), add.names = c("NCI_2","AD_2"))
cellchat_merged3 <- mergeCellChat(list(NCI_niche3, AD_niche3), add.names = c("NCI_3","AD_3"))
saveRDS(cellchat_merged, file = "cellchat_merged_full.rds")
saveRDS(cellchat_merged1, file = "cellchat_merged_k3_n1.rds")
saveRDS(cellchat_merged2, file = "cellchat_merged_k3_n2.rds")
saveRDS(cellchat_merged3, file = "cellchat_merged_k3_n3.rds")

# merge and save normalized objects
cellchat_merged_norm <- mergeCellChat(list(cellchat_nci_norm, cellchat_ad_norm), add.names = c("norm_NCI","norm_AD"))
cellchat_merged_norm1 <- mergeCellChat(list(cellchat_nci_niche1_norm, cellchat_ad_niche1_norm), add.names = c("norm_NCI_1","norm_AD_1"))
cellchat_merged_norm2 <- mergeCellChat(list(cellchat_nci_niche2_norm, cellchat_ad_niche2_norm), add.names = c("norm_NCI_2","norm_AD_2"))
cellchat_merged_norm3 <- mergeCellChat(list(cellchat_nci_niche3_norm, cellchat_ad_niche3_norm), add.names = c("norm_NCI_3","norm_AD_3"))
saveRDS(cellchat_merged_norm, file = "cellchat_merged_norm_full.rds")
saveRDS(cellchat_merged_norm1, file = "cellchat_merged_norm_k3_n1.rds")
saveRDS(cellchat_merged_norm2, file = "cellchat_merged_norm_k3_n2.rds")
saveRDS(cellchat_merged_norm3, file = "cellchat_merged_norm_k3_n3.rds")

# heatmaps by count (unnormalized) or weight (normalized; Extended Data Figure 6)
netVisual_heatmap(cellchat_merged, measure = "count", color.use = major)
netVisual_heatmap(cellchat_merged_norm, measure = "weight", color.use = major) 
netVisual_heatmap(cellchat_merged1, measure = "count", color.use = major) 
netVisual_heatmap(cellchat_merged_norm1, measure = "weight", color.use = major) 
netVisual_heatmap(cellchat_merged2, measure = "count", color.use = major) 
netVisual_heatmap(cellchat_merged_norm2, measure = "weight", color.use = major) 
netVisual_heatmap(cellchat_merged3, measure = "count", color.use = major) 
netVisual_heatmap(cellchat_merged_norm3, measure = "weight", color.use = major) 

### AD-differential heatmap (Figure 6i)

# 1. Define cellchat objects and AD/NCI pairs
cellchat_list <- list(
  NCI     = cellchat_nci_norm,
  AD      = cellchat_ad_norm,
  NCI_n1  = cellchat_nci_niche1_norm,
  AD_n1   = cellchat_ad_niche1_norm,
  NCI_n2  = cellchat_nci_niche2_norm,
  AD_n2   = cellchat_ad_niche2_norm,
  NCI_n3  = cellchat_nci_niche3_norm,
  AD_n3   = cellchat_ad_niche3_norm
)

condition_pairs <- list(
  Full   = c("AD", "NCI"),
  Niche1 = c("AD_n1", "NCI_n1"),
  Niche2 = c("AD_n2", "NCI_n2"),
  Niche3 = c("AD_n3", "NCI_n3")
)

# 2. Get all pathway names
all_pathways <- unique(unlist(lapply(cellchat_list, function(cellchat) {
  if (!is.null(cellchat@netP$prob)) {
    return(dimnames(cellchat@netP$prob)[[3]])
  } else {
    return(NULL)
  }
})))

# 3. Compute mean netP probabilities per condition for each pathway
pathway_probs_df <- data.frame(row.names = all_pathways)

for (name in names(cellchat_list)) {
  cellchat <- cellchat_list[[name]]
  if (!is.null(cellchat@netP$prob)) {
    prob_array <- cellchat@netP$prob
    pathway_probs <- apply(prob_array, 3, function(x) mean(x, na.rm = TRUE))
    pathway_probs_filtered <- pathway_probs[pathway_probs >= 0]  # Filter near-zero
    matched_probs <- rep(NA, length(all_pathways))
    names(matched_probs) <- all_pathways
    matched_probs[names(pathway_probs_filtered)] <- pathway_probs_filtered
    pathway_probs_df[[name]] <- matched_probs
  } else {
    pathway_probs_df[[name]] <- rep(NA, length(all_pathways))
  }
}

# 4. Compute difference matrix (AD - NCI)
diff_df <- data.frame(row.names = all_pathways)

for (pair_name in names(condition_pairs)) {
  ad_col  <- condition_pairs[[pair_name]][1]
  nci_col <- condition_pairs[[pair_name]][2]
  diff <- pathway_probs_df[[ad_col]] - pathway_probs_df[[nci_col]]
  diff_df[[pair_name]] <- diff
}

# 5. Select top 20 pathways based on absolute difference in the full comparison
diff_df$diff_full <- diff_df$Full
diff_df$abs_diff_full <- abs(diff_df$diff_full)
top20_pathways <- rownames(diff_df)[order(-diff_df$abs_diff_full)][1:20]
# or manually define
sig_subset <- rev(c("COLLAGEN",'APP',"LAMININ",'FN1','FGF','THBS','NOTCH','ADM','MIF','VEGF','ADGRG','PDGF','CXCL','IGF','CSF','VCAM','RA','Glutamate','SPP1','ApoE','WNT'))
pathways <- sig_subset # or manually define

# 6. Sort top 20 pathways by AD-NCI difference in full dataset
pathways_sorted <- pathways[order(-diff_df[pathways, "diff_full"])]

# 7. Prepare matrix for heatmap
diff_mat <- as.matrix(diff_df[pathways_sorted, names(condition_pairs)])

# 8. Plot heatmap (Fig. 6k)
max_val <- max(abs(diff_mat), na.rm = TRUE)/3

pheatmap(
  diff_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Top Pathways - Difference (AD - NCI)",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-max_val, max_val, length.out = 51),
  fontsize_row = 10
)


# load seurat object
ChP_filt <- readRDS("data_objects/ChP_filt_082025.rds") # small object, without images or SCT assay

# cell type colors
major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B") 

# merged norm CC objects; full or niche1-3
cellchat_merged_norm <- readRDS(file = "data_objects/cellchat_merged_norm_full.rds")
cellchat_merged_norm1 <- readRDS(file = "data_objects/cellchat_merged_norm_k3_n1.rds") # niche1
cellchat_merged_norm2 <- readRDS(file = "data_objects/cellchat_merged_norm_k3_n2.rds") # niche2
cellchat_merged_norm3 <- readRDS(file = "data_objects/cellchat_merged_norm_k3_n3.rds") # niche3

# generate 'net' dfs based on CC
Ex_ChP_merge_CC <- identifyOverExpressedGenes(cellchat_merged_norm, group.dataset = "datasets", pos.dataset = "norm_AD", features.name = "norm_AD", only.pos = FALSE)
net <- netMappingDEG(Ex_ChP_merge_CC, features.name = "norm_AD")
Ex_ChP_merge_CC1 <- identifyOverExpressedGenes(cellchat_merged_norm1, group.dataset = "datasets", pos.dataset = "norm_AD_1", features.name = "norm_AD_1", only.pos = FALSE)
net1 <- netMappingDEG(Ex_ChP_merge_CC1, features.name = "norm_AD_1")
Ex_ChP_merge_CC2 <- identifyOverExpressedGenes(cellchat_merged_norm3, group.dataset = "datasets", pos.dataset = "norm_AD_3", features.name = "norm_AD_3", only.pos = FALSE)
net2 <- netMappingDEG(Ex_ChP_merge_CC2, features.name = "norm_AD_3")
Ex_ChP_merge_CC3 <- identifyOverExpressedGenes(cellchat_merged_norm2, group.dataset = "datasets", pos.dataset = "norm_AD_2", features.name = "norm_AD_2", only.pos = FALSE)
net3 <- netMappingDEG(Ex_ChP_merge_CC3, features.name = "norm_AD_2")

# 1. subset seurat object by niche
ChP_filt_nano_1 <- subset(ChP_filt_nano, subset = niche_k3 == "1")
ChP_filt_nano_2 <- subset(ChP_filt_nano, subset = niche_k3 == "2")
ChP_filt_nano_3 <- subset(ChP_filt_nano, subset = niche_k3 == "3")

# add dimnames back to seurat objects; subsetting can cause them to become null
# 2. Get the gene and cell names
genes <- rownames(ChP_filt_nano)
cells1 <- colnames(ChP_filt_nano_1)
cells2 <- colnames(ChP_filt_nano_2)
cells3 <- colnames(ChP_filt_nano_3)

# 3. Assign gene and cell names to the counts layer
ChP_filt_nano_1@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells1)
ChP_filt_nano_2@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells2)
ChP_filt_nano_3@assays[["Nanostring"]]@layers[["counts"]]@Dimnames <- list(genes, cells3)

# 4. ensure AD status is stored as '1' or '2' and as a character
ChP_filt_nano$AD_status <- recode(ChP_filt_nano$AD_status, "AD" = "1", "NCI" = "2")
ChP_filt_nano_1$AD_status <- recode(ChP_filt_nano_1$AD_status, "AD" = "1", "NCI" = "2")
ChP_filt_nano_2$AD_status <- recode(ChP_filt_nano_2$AD_status, "AD" = "1", "NCI" = "2")
ChP_filt_nano_3$AD_status <- recode(ChP_filt_nano_3$AD_status, "AD" = "1", "NCI" = "2")

ChP_filt_nano$AD_status <- as.character(ChP_filt_nano$AD_status)
ChP_filt_nano_1$AD_status <- as.character(ChP_filt_nano_1$AD_status)
ChP_filt_nano_2$AD_status <- as.character(ChP_filt_nano_2$AD_status)
ChP_filt_nano_3$AD_status <- as.character(ChP_filt_nano_3$AD_status)


# Run for each niche, or full dataset
n1_results <- process_niche(ChP_filt_nano_1, net1, niche_name = "1")
n2_results <- process_niche(ChP_filt_nano_2, net2, niche_name = "2")
n3_results <- process_niche(ChP_filt_nano_3, net3, niche_name = "3")
nAll_results <- process_niche(ChP_filt_nano, net, niche_name = "all")

# save
saveRDS(nAll_results, file = "data_objects/Full_CC_DEG_net.rds")
saveRDS(n1_results, file = "data_objects/n1_CC_DEG_net.rds")
saveRDS(n2_results, file = "data_objects/n2_CC_DEG_net.rds")
saveRDS(n3_results, file = "data_objects/n3_CC_DEG_net.rds")

# start from here after having saved updated net results
nAll_results <- readRDS(file = "data_objects/Full_CC_DEG_net.rds")
n1_results <- readRDS(file = "data_objects/n1_CC_DEG_net.rds")
n2_results <- readRDS(file = "data_objects/n2_CC_DEG_net.rds")
n3_results <- readRDS(file = "data_objects/n3_CC_DEG_net.rds")

# load CC objects
cellchat_nci_norm <- readRDS(file = "data_objects/cellchat_nci_norm.rds")
cellchat_ad_norm <- readRDS(file = "data_objects/cellchat_ad_norm.rds")
cellchat_nci_niche1_norm <- readRDS(file = "data_objects/cellchat_nci_k3_niche1_norm.rds")
cellchat_ad_niche1_norm <- readRDS(file = "data_objects/cellchat_ad_k3_niche1_norm.rds")
cellchat_nci_niche2_norm <- readRDS(file = "data_objects/cellchat_nci_k3_niche2_norm.rds")
cellchat_ad_niche2_norm <- readRDS(file = "data_objects/cellchat_ad_k3_niche2_norm.rds")
cellchat_nci_niche3_norm <- readRDS(file = "data_objects/cellchat_nci_k3_niche3_norm.rds")
cellchat_ad_niche3_norm <- readRDS(file = "data_objects/cellchat_ad_k3_niche3_norm.rds")

ggsave("figures/full_combined_logFC_plot.pdf", plot = nAll_results$plot_combined, width = 6, height = 5)
ggsave("figures/niche1_combined_logFC_plot.pdf", plot = n1_results$plot_combined, width = 6, height = 5)
ggsave("figures/niche2_combined_logFC_plot.pdf", plot = n2_results$plot_combined, width = 6, height = 5)
ggsave("figures/niche3_combined_logFC_plot.pdf", plot = n3_results$plot_combined, width = 6, height = 5)



# Plot the distribution of % cells expressing ligand or receptor (overlaid)
plot_pct_expressing_distribution(nAll_results$updated_net)
plot_pct_expressing_distribution(n1_results$updated_net)
plot_pct_expressing_distribution(n2_results$updated_net)
plot_pct_expressing_distribution(n3_results$updated_net)

# generate significance matrices to be used when plotting AD-diff centrality heatmaps
# STEP 1: Filter net for valid logFC and pathway info
net_filtered <- nAll_results$updated_net %>%
  filter(
    !is.na(combined_logFC),
    !is.na(pathway_name),
    prob > 5e-07, # optionally, filter out low prob interactions; this is normalized value
    abs(combined_logFC) > 0.05, # optionally, filter out noisy/non-directional effects
    ligand.pct_expressing > 5, # expressed in >5% of cells (ligand-sender) 
    receptor.pct_expressing > 5) # expressed in >5% of cells (receptor-target) 

# prob and logFC filters - from 13K to 10K
# ligand/receptor exp filters - only an additional 34 interactions

# STEP 2: Define cell types and pathways of interest
cell_types_of_interest <- c("Epithelial", "Fibroblast", "Endothelial", "Mural", "Immune")
pathways <- unique(net_filtered$pathway_name)

# STEP 3A: Run signed-rank tests pathway × celltype (source/target collapsed)
results_list <- list()
for (pw in pathways) {
  for (ct in cell_types_of_interest) {
    df_sub <- net_filtered %>%
      filter(pathway_name == pw, source == ct | target == ct)
    
    if (nrow(df_sub) == 0) next
    
    n_pairs <- nrow(df_sub)
    median_fc <- median(df_sub$combined_logFC, na.rm = TRUE)
    
    # Wilcoxon test for deviation from 0
    p_val <- if (n_pairs >= 3) {
      wilcox.test(df_sub$combined_logFC, mu = 0, alternative = "two.sided")$p.value
    } else {
      NA
    }
    
    results_list[[length(results_list) + 1]] <- data.frame(
      pathway_name = pw,
      cell_type = ct,
      n_pairs = n_pairs,
      median_logFC = median_fc,
      p_signed_rank = p_val
    )
  }
}
combined_stats <- bind_rows(results_list)
# STEP 3B: Adjust p-values for multiple testing
combined_stats$FDR <- p.adjust(combined_stats$p_signed_rank, method = "fdr")
write.csv(combined_stats, file="data_csv/full_GSEA_pathway_cell-type_stats.csv")

# STEP 4: Run signed-rank tests pathway x sender x receiver
results_list <- list()
for (pw in unique(net_filtered$pathway_name)) {
  df_pw <- net_filtered %>% filter(pathway_name == pw)
  
  # Get all unique source–target pairs for this pathway
  pairs <- unique(df_pw[, c("source", "target")])
  
  for (i in seq_len(nrow(pairs))) {
    src <- pairs$source[i]
    tgt <- pairs$target[i]
    
    df_sub <- df_pw %>% filter(source == src, target == tgt)
    
    if (nrow(df_sub) == 0) next
    
    n_pairs <- nrow(df_sub)
    median_fc_lig <- median(df_sub$ligand.logFC, na.rm = TRUE)
    median_fc_rec <- median(df_sub$receptor.logFC, na.rm = TRUE)
    median_fc <- median(df_sub$combined_logFC, na.rm = TRUE)
    
    p_val <- if (n_pairs >= 3) {
      wilcox.test(df_sub$combined_logFC, mu = 0)$p.value
    } else {
      NA
    }
    
    results_list[[length(results_list) + 1]] <- data.frame(
      source = src,
      target = tgt,
      pathway_name = pw,
      n_pairs = n_pairs,
      median_ligand_logFC = median_fc_lig,
      median_receptor_logFC = median_fc_rec,
      median_logFC = median_fc,
      p_signed_rank = p_val
    )
  }
}
combined_stats_full <- bind_rows(results_list)
combined_stats_full$FDR <- p.adjust(combined_stats_full$p_signed_rank, method = "fdr")
# above version is for supplemental netP table

# STEP 5: Filter and prepare for plotting
sig_subset <- rev(c("COLLAGEN",'APP',"LAMININ",'FN1','FGF','THBS','NOTCH','ADM','MIF','VEGF','ADGRG','PDGF','CXCL','IGF','CSF','VCAM','RA','Glutamate','SPP1','ApoE','WNT'))

filtered_stats <- combined_stats %>%
  filter(pathway_name %in% sig_subset) %>%
  mutate(pathway_name = factor(pathway_name, levels = sig_subset)) %>%
  filter(FDR < 0.05 & abs(median_logFC) > 0.15)

# plot summary data
ggplot(combined_stats, aes(x = median_logFC, y = -log10(FDR), color = cell_type)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  scale_color_manual(values = major) +
  theme_minimal() +
  labs(
    x = "Median Combined logFC",
    y = "-log10(FDR)",
    color = "Cell Type",
    title = "Full: Pathway-level Summary (per Cell Type)"
  ) +
  theme(
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

# STEP 6: Create full grid for annotation
full_df <- expand.grid(
  pathway_name = factor(sig_subset, levels = sig_subset),
  cell_type = cell_types_of_interest,
  stringsAsFactors = FALSE
) %>%
  left_join(filtered_stats %>% select(pathway_name, cell_type, FDR), by = c("pathway_name", "cell_type")) %>%
  mutate(
    FDR = ifelse(is.na(FDR), 1, FDR),
    signif = case_when(
      FDR < 0.01  ~ "*",
      TRUE        ~ ""
    )
  )

# Convert to matrix format: pathways (rows) x cell types (columns)
full_df <- full_df %>%
  filter(!is.na(pathway_name)) 

signif_mat <- full_df %>%
  select(pathway_name, cell_type, signif) %>%
  pivot_wider(names_from = cell_type, values_from = signif, values_fill = "") %>%
  column_to_rownames("pathway_name") %>%
  as.matrix()
signif_mat_use <- t(signif_mat)

# STEP 7: Plot heatmap of asterisks
ggplot(full_df, aes(x = factor(cell_type, levels = cell_types_of_interest),
                    y = pathway_name)) +
  geom_text(aes(label = signif), size = 6) +
  theme_minimal() +
  labs(
    x = "Cell Type",
    y = "Pathway",
    title = "Full: Pathway-Level Cell-Type Significance (Wilcoxon, FDR-adjusted)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## custom heatmap with asterisks
# prep matrix (load function first)
mat1 <- prep_mat(object=cellchat_nci_norm, slot.name="netP", thresh=0.05, pattern="all")
mat3 <- prep_mat(object=cellchat_ad_norm, slot.name="netP", thresh=0.05, pattern="all")
mat <- bind_rows(mat1, mat3)

curr_dir <- '~'
netAnalysis_signalingRole_heatmap_diff_or_ratio_comp(
  mat,
  pattern = "all",
  signaling = sig_subset,
  filename = "AD_vs_NCI_full_0.15_0.01",
  comparison.type = "diff",
  color.use = c(Epithelial = "#C48F45", Fibroblast = "#3249A6", Endothelial = "#926392", Mural = "#261132", Immune = "#26532B"),
  significance.matrix = signif_mat_use
)

# incorporate differential centrality into asterisks as well (direction must match with average logFC)
# Ensure common pathways and cell types
common_pathways <- intersect(colnames(mat3), colnames(mat1))
common_celltypes <- intersect(rownames(mat3), rownames(mat1))

# Subset to shared content
mat_ad_sub <- mat3[common_celltypes, common_pathways]
mat_nci_sub <- mat1[common_celltypes, common_pathways]

mat_ad_sub <- mat_ad_sub[, !(colnames(mat_ad_sub) %in% c("group"))]
mat_nci_sub <- mat_nci_sub[, !(colnames(mat_nci_sub) %in% c("group"))]

# Compute difference: AD − NCI
centrality_diff <- mat_ad_sub - mat_nci_sub
rownames(centrality_diff) <- cell_types_of_interest

# filter for subset of pathways of interest
centrality_diff_df <- as.data.frame(centrality_diff) %>%
  rownames_to_column("cell_type") %>%
  pivot_longer(-cell_type, names_to = "pathway_name", values_to = "centrality_diff") %>%
  filter(pathway_name %in% sig_subset)

# revised significance, incorporating both GSEA and centrality diff
full_df <- expand.grid(
  pathway_name = factor(sig_subset, levels = sig_subset),
  cell_type = cell_types_of_interest,
  stringsAsFactors = FALSE
) %>%
  left_join(filtered_stats %>% select(pathway_name, cell_type, FDR, median_logFC), 
            by = c("pathway_name", "cell_type")) %>%
  left_join(centrality_diff_df %>% select(pathway_name, cell_type, centrality_diff), 
            by = c("pathway_name", "cell_type")) %>%
  mutate(
    centrality_sig = abs(centrality_diff) > 0.0001,       # your chosen threshold
    deg_sig = FDR < 0.05,
    direction_match = sign(median_logFC) == sign(centrality_diff),
    show_star = centrality_sig & deg_sig & direction_match,
    signif = case_when(
      show_star & FDR < 0.01   ~ "*",
      TRUE                    ~ ""
    )) %>%
  mutate(
    pathway_name = factor(pathway_name, levels = sig_subset)  # reset order explicitly
  )

full_df <- full_df %>%
  filter(!is.na(pathway_name)) 

# STEP 8: Plot heatmap of asterisks
ggplot(full_df, aes(x = factor(cell_type, levels = cell_types_of_interest),
                    y = pathway_name)) +
  geom_text(aes(label = signif), size = 6) +
  theme_minimal() +
  labs(
    x = "Cell Type",
    y = "Pathway",
    title = "Full: Pathway-Level Cell-Type Significance (Wilcoxon, FDR-adjusted), v2"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Convert to matrix format: pathways (rows) x cell types (columns)
signif_mat <- full_df %>%
  select(pathway_name, cell_type, signif) %>%
  pivot_wider(names_from = cell_type, values_from = signif, values_fill = "") %>%
  column_to_rownames("pathway_name") %>%
  as.matrix()
signif_mat_use <- t(signif_mat)

# optionally, re-generate with revised sig (or just copy asterisks from asterisk plot above)
# Figure 6j
netAnalysis_signalingRole_heatmap_diff_or_ratio_comp(
  mat,
  pattern = "all",
  signaling = sig_subset,
  filename = "AD_vs_NCI_full_0.2_0.01_v3",
  comparison.type = "diff",
  color.use = c(Epithelial = "#C48F45", Fibroblast = "#3249A6", Endothelial = "#926392", Mural = "#261132", Immune = "#26532B"),
  significance.matrix = signif_mat_use
)

# repeat steps 1-8 above for niches 1-3

### prob vs. logFC scatter (extra supplemental plot (not included in manuscript))

# Subset for a specific pathway, e.g., COLLAGEN
pathway_of_interest <- "COLLAGEN"

# convert net for full/niche to a df
net_filtered <- nAll_results$updated_net
net1_filtered <- n1_results$updated_net
net2_filtered <- n2_results$updated_net
net3_filtered <- n3_results$updated_net

# Step 1: Reshape to get AD and NCI side-by-side
net_diff <- net_filtered %>%
  pivot_wider(
    names_from = datasets,
    values_from = c(prob, combined_logFC),
    names_sep = "."
  ) %>%
  mutate(
    delta_prob = prob.norm_AD - prob.norm_NCI,
    combined_logFC = rowMeans(
      cbind(combined_logFC.norm_AD, combined_logFC.norm_NCI),
      na.rm = TRUE
    ),
    relative_delta = delta_prob / max(abs(delta_prob), na.rm = TRUE)
  )

net1_diff <- net1_filtered %>%
  pivot_wider(
    names_from = datasets,
    values_from = c(prob, combined_logFC),
    names_sep = "."
  ) %>%
  mutate(
    delta_prob = prob.norm_AD_1 - prob.norm_NCI_1,
    combined_logFC = rowMeans(
      cbind(combined_logFC.norm_AD_1, combined_logFC.norm_NCI_1),
      na.rm = TRUE
    ),
    relative_delta = delta_prob / max(abs(delta_prob), na.rm = TRUE)
  )

net2_diff <- net2_filtered %>%
  pivot_wider(
    names_from = datasets,
    values_from = c(prob, combined_logFC),
    names_sep = "."
  ) %>%
  mutate(
    delta_prob = prob.norm_AD_3 - prob.norm_NCI_3,
    combined_logFC = rowMeans(
      cbind(combined_logFC.norm_AD_3, combined_logFC.norm_NCI_3),
      na.rm = TRUE
    ),
    relative_delta = delta_prob / max(abs(delta_prob), na.rm = TRUE)
  )

net3_diff <- net3_filtered %>%
  pivot_wider(
    names_from = datasets,
    values_from = c(prob, combined_logFC),
    names_sep = "."
  ) %>%
  mutate(
    delta_prob = prob.norm_AD_2 - prob.norm_NCI_2,
    combined_logFC = rowMeans(
      cbind(combined_logFC.norm_AD_2, combined_logFC.norm_NCI_2),
      na.rm = TRUE
    ),
    relative_delta = delta_prob / max(abs(delta_prob), na.rm = TRUE)
  )

# filter and apply thresholds for highlighting
plot_df <- net_filtered %>%
  filter(pathway_name == pathway_of_interest) %>%
  mutate(highlight = ifelse(prob > 2.5e-5 & abs(combined_logFC) > 0.2, "High Activity", "Other"))

ggplot(plot_df, aes(x = combined_logFC, y = prob)) +
  geom_point(aes(color = highlight), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("High Activity" = "red", "Other" = "gray70")) +
  geom_hline(yintercept = 2.5e-5, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray50") +
  theme_minimal() +
  labs(
    title = paste0("Full - Pathway Activity: ", pathway_of_interest),
    subtitle = "Highlighting high-probability, high-logFC interactions",
    x = "Combined logFC (ligand + receptor)",
    y = "Interaction Probability",
    color = "Interaction Type"
  )



pathways_to_plot <- c("COLLAGEN", "APP", "LAMININ", "FN1")

plot_pathways_pdf(
  pathways = pathways_to_plot,
  major = major,
  label_LR = TRUE,
  output_file = "figures/LR_interaction_grid_top5_070925.pdf"
)


plot_rel_pathways_pdf(
  pathways = pathways_to_plot3,
  major = major,
  label_LR = TRUE,
  output_file = "figures/LR_interaction_grid_relative_0.25_070925_3.pdf"
)


