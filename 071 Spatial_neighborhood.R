###
### 071 Spatial_neighborhood.R
### This script was written by Dr. Denis R Avey
# Purpose: Compute and cluster spatial neighborhoods based on neighbor-averaged expression.
# Dependencies:
source("/070 Spatial_neighborhood_util", verbose = FALSE)
curr_dir <- "/ChP_ST" #output results to specific folder


# load ChP_filt data and parse relevant data
ChP_filt <- readRDS("ChP_filt.rds")

annot_meta <- ChP_filt@meta.data
rownames(annot_meta) <- names(ChP_filt@active.ident) # non-unique values detected in @meta.data, so use rownames of active.ident instead
counts <- ChP_filt@assays[["Nanostring"]]@layers[["counts"]] # raw counts
rownames(counts) <- rownames(ChP_filt@assays[["SCT_full"]]@counts) 
colnames(counts) <- colnames(ChP_filt@assays[["SCT_full"]]@counts)
celltype <- as.factor(annot_meta$newMajorcelltype)

xy1 <- ChP_filt@images[["Slide.1"]]@boundaries[["centroids"]]@coords
xy2 <- ChP_filt@images[["Slide.2"]]@boundaries[["centroids"]]@coords
xy <- rbind(xy1,xy2)
rownames(xy) <- annot_meta$unique_cellID
rm(xy1)
rm(xy2)
saveRDS(xy, file="ChP_xy.rds")
saveRDS(ChP_filt, file= "ChP_filt.rds")
rm(ChP_filt) # remove to save RAM and reload later
gc()

# define neighbors using a K-nearest approach:
neighbors.nearest25 <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 25) 
# define using a radius-based approach:
# neighbors.radiusbased <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 500)

rm(xy)

# the output is a sparse matrix of cells * cells:

# compare the number of neighbors found by each approach:

# In “Mitigating autocorrelation during spatially resolved transcriptomics data analysis”, Maher et al. (2023) describe an inconvenient tendency of spatial context matrices: because neighboring cells have largely the same neighbors, their entries in the spatial context matrix are correlated. This correlation between neighbors proves a substantial barrier to distance-based analyses like UMAP or Leiden clustering, producing UMAPs where all points fall in a highly-connected blob and generally poor Leiden performance. (However, for most analyses, correlation between neighboring cells’ spatial context vectors doesn’t seem to have much impact.) They propose that by defining each cell’s neighborhood as a random subset of its nearest neighbors, they can largely break this correlation between neighbors. They released a python toolkit for this, SPIN.

# For R coders, here’s how you would get a neighborhood matrix with random subsetting:
subsetted_neighbors10 <- subsampleNeighborsByRow(neighbors = neighbors.nearest25, p = 0.4) # change p to change % downsampled, this downsamples to a random 10-cell subset of the 25 nearest neighbors

# optionally, include the cell itself in the neighborhood; this makes sense for downstream Cellchat, especially if we want to interrogate cell-cell interactions among rare cell types
diag(neighbors.nearest25) <- 1  # Add the cell itself to the diagonal
diag(subsetted_neighbors10) <- 1  # Add the cell itself to the diagonal


# re-format to get spatial context matrix
counts <- t(counts) #transpose so that cells are rows and genes are columns
rownames(neighbors.nearest25) <- rownames(counts) # add rownames to neighbors.nearest25
rownames(subsetted_neighbors10) <- rownames(counts) # add rownames to subsetted_neighbors

# mean neighborhood expression:
spatialcontext_expression <- get_neighborhood_expression(counts, neighbors.nearest25) # highly RAM-intensive
spatialcontext_expression_sub10 <- get_neighborhood_expression(counts, subsetted_neighbors10)
gc()
# Convert the matrix to a Seurat Assay
ChP_filt <- readRDS("ChP_filt_032825.rds")
cell_neighborhoods_25 <- CreateAssayObject(data = t(spatialcontext_expression))
ChP_filt[["SpatialNeighborhood"]] <- cell_neighborhoods_25

# Repeat for subsetted neighbors
cell_neighborhoods_25_sub10 <- CreateAssayObject(data = t(spatialcontext_expression_sub10))
ChP_filt[["SpatialNeighborhood_sub10"]] <- cell_neighborhoods_25_sub10

# remove spatial matrices to save RAM
rm(spatialcontext_expression)
rm(spatialcontext_expression_sub10)
gc()

# Copy the data slot to the counts slot
ChP_filt[["SpatialNeighborhood"]]@counts <- ChP_filt[["SpatialNeighborhood"]]@data
ChP_filt[["SpatialNeighborhood_sub10"]]@counts <- ChP_filt[["SpatialNeighborhood_sub10"]]@data
# Check counts slot
# Check data slot

# Convert 'data' slot to sparse matrix
ChP_filt[["SpatialNeighborhood"]]@data <- as(as.matrix(ChP_filt[["SpatialNeighborhood"]]@data), "dgCMatrix")
# Convert 'counts' slot to sparse matrix
ChP_filt[["SpatialNeighborhood"]]@counts <- as(as.matrix(ChP_filt[["SpatialNeighborhood"]]@counts), "dgCMatrix")

# Convert 'data' slot to sparse matrix
ChP_filt[["SpatialNeighborhood_sub10"]]@data <- as(as.matrix(ChP_filt[["SpatialNeighborhood_sub10"]]@data), "dgCMatrix")
ChP_filt[["SpatialNeighborhood_sub10"]]@counts <- as(as.matrix(ChP_filt[["SpatialNeighborhood_sub10"]]@counts), "dgCMatrix")

# Now that the SpatialNeighborhood assay is part of the Seurat object, you can perform standard Seurat workflows on it, such as scaling, clustering, or visualization.

# Normalize and Scale Data
ChP_filt <- NormalizeData(ChP_filt, assay = "SpatialNeighborhood")
ChP_filt <- ScaleData(ChP_filt, assay = "SpatialNeighborhood")
# Find Variable Features
ChP_filt <- FindVariableFeatures(ChP_filt, assay = "SpatialNeighborhood")
# Run PCA or Other Reductions
ChP_filt <- RunPCA(ChP_filt, assay = "SpatialNeighborhood", reduction.name = "spatial_pca_25_NEW")

# Run Harmony batch correction on PCA
ChP_filt <- RunHarmony(
  ChP_filt,
  group.by.vars = c("Slide", "projID","AD_status"),  # Variables indicating slide number and projID; test including AD here as well
  reduction.use = "spatial_pca_25_NEW",  # Use the correct spatial PCA reduction
  dims.use = 1:50,  # Adjust based on PCA exploration
  reduction.save = "spatial_harmony_25_NEW"  # Specify the name for the Harmony reduction
)

# Re-run UMAP on Harmony embeddings
ChP_filt <- RunUMAP(ChP_filt, reduction = "spatial_harmony_25_NEW", dims = 1:50, reduction.name = "spatial_harmony_25_umap_NEW")

# repeat for subsampled neighborhood matrix, 10-cells
# Normalize and Scale Data
ChP_filt <- NormalizeData(ChP_filt, assay = "SpatialNeighborhood_sub10")
ChP_filt <- ScaleData(ChP_filt, assay = "SpatialNeighborhood_sub10")
# Find Variable Features
ChP_filt <- FindVariableFeatures(ChP_filt, assay = "SpatialNeighborhood_sub10")
# Run PCA or Other Reductions
ChP_filt <- RunPCA(ChP_filt, assay = "SpatialNeighborhood_sub10", reduction.name = "spatial_pca_25_sub10")

# Run Harmony batch correction on PCA
ChP_filt <- RunHarmony(
  ChP_filt,
  group.by.vars = c("Slide", "projID","AD_status"),  # Variables indicating slide number and projID; test including AD here as well
  reduction.use = "spatial_pca_25_sub10",  # Use the correct spatial PCA reduction
  dims.use = 1:50,  # Adjust based on PCA exploration
  reduction.save = "spatial_harmony_25_sub10"  # Specify the name for the Harmony reduction
)

# Re-run UMAP on Harmony embeddings
ChP_filt <- RunUMAP(ChP_filt, reduction = "spatial_harmony_25_sub10", dims = 1:50, reduction.name = "spatial_harmony_25_umap_sub10")

saveRDS(ChP_filt, file = "ChP_filt_neighbors.rds") # save new Seurat object containing spatial expression matrix/pca/harmony/umap with low-quality and Parenchyma filtered out

# kmeans clustering (using 'sub10' version because there is less spatial autocorrelation between nearby cells)
ChP_filt_sub10 <- readRDS("ChP_filt_neighbors.rds")

# Extract Harmony embeddings
harmony_embeddings_25_sub10 <- ChP_filt_sub10@reductions[["spatial_harmony_25_sub10"]]@cell.embeddings

# justify k number
wcss <- sapply(2:12, function(k) {
  kmeans(harmony_embeddings_25_sub10, centers = k, nstart = 25)$tot.withinss
})

elbow_plot <- data.frame(k = 2:12, wcss = wcss)

ggplot(elbow_plot, aes(x = k, y = wcss)) +
  geom_line() + geom_point() +
  labs(title = "Elbow Method for Optimal k", x = "Number of Clusters (k)", y = "WCSS")

# Run K-means clustering (test multiple versions)
set.seed(123)
kmeans_result2_sub10 <- kmeans(harmony_embeddings_25_sub10, centers = 2)
ChP_filt_sub10$niche_k2 <- factor(kmeans_result2_sub10$cluster)

kmeans_result3_sub10 <- kmeans(harmony_embeddings_25_sub10, centers = 3)
ChP_filt_sub10$niche_k3 <- factor(kmeans_result3_sub10$cluster)

kmeans_result4_sub10 <- kmeans(harmony_embeddings_25_sub10, centers = 4)
ChP_filt_sub10$niche_k4 <- factor(kmeans_result4_sub10$cluster)

# define niche colors
niche_col <- c(
  "1"    = "#33a02c",   
  "2"    = "#1f78b4",
  "3"    = "#ffb000")

# PCA plot of neighborhood clusters (Figure 6d: left)
DimPlot(ChP_filt_sub10, reduction = "spatial_pca_25_sub10", cols = niche_col, group.by = "niche_k3")

# Plot niche clusters in the FOV
DefaultAssay(ChP_filt_sub10) <- "SpatialNeighborhood_sub10" # Ensure the correct assay is set
DefaultBoundary(ChP_filt_sub10@images$Slide.1.89) <- 'segmentation'

# ImageDimPlot of FOV89, colored by niche (Figure 6d: right)
ImageDimPlot(ChP_filt_sub10, fov = "Slide.1.89", cols=niche_col, alpha = 1, group.by = "niche_k3", dark.background = FALSE, border.color = "black",coord.fixed = FALSE, border.size = 0.2)

# major cell type colors
major <- c("Epithelial"="#C48F45", "Fibroblast"="#3249A6", "Endothelial"="#926392", "Mural"="#261132", "Immune"="#26532B")

# proportion plots, grouped by niche
pdf("proportion_plots_sub10_clusters.pdf")
celltype_props <- ChP_filt_sub10@meta.data %>%
  group_by(niche_k3, newMajorcelltype) %>%
  summarise(count = n()) %>%
  group_by(niche_k3) %>%
  mutate(proportion = count / sum(count))

# Plot as stacked barplot
ggplot(celltype_props, aes(x = niche_k3, y = proportion, fill = newMajorcelltype)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +  # Add black border around bars
  theme_minimal() +
  labs(title = "Cell-Type Proportions by Spatial Cluster Kmeans-3", x = "Spatial Cluster", y = "Proportion") +
  scale_fill_manual(values = major)

# AD_status, gender, batch comparison #

# Summarize cluster proportions by AD_status
cluster_summary <- ChP_filt_sub10@meta.data %>%
  group_by(AD_status, niche_k3) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(AD_status) %>%
  mutate(proportion = count / sum(count))

# Plot cluster proportions by AD status
ggplot(cluster_summary, aes(x = niche_k3, y = proportion, fill = AD_status)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "Spatial25 Cluster Proportions by AD Status", x = "Spatial Cluster", y = "Proportion") +
  scale_fill_manual(values = c("NCI" = "lightblue", "AD" = "salmon"))

# Repeat for gender 
cluster_gender <- ChP_filt_sub10@meta.data %>%
  group_by(gender, niche_k3) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(gender) %>%
  mutate(proportion = count / sum(count))

ggplot(cluster_gender, aes(x = niche_k3, y = proportion, fill = gender)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "Spatial25 Cluster Proportions by Gender", x = "Spatial Cluster", y = "Proportion")

# Repeat for Slide (technical batch) 
cluster_slide <- ChP_filt_sub10@meta.data %>%
  group_by(Slide, niche_k3) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Slide) %>%
  mutate(proportion = count / sum(count))

ggplot(cluster_slide, aes(x = niche_k3, y = proportion, fill = Slide)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "Spatial25 Cluster Proportions by Slide", x = "Spatial Cluster", y = "Proportion")

# Repeat for individual (arb_ID) 
cluster_arb_ID <- ChP_filt_sub10@meta.data %>%
  group_by(arb_ID, niche_k3) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(arb_ID) %>%
  mutate(proportion = count / sum(count))

# Plot the proportions by arb_ID
ggplot(cluster_arb_ID, aes(x = niche_k3, y = proportion, fill = arb_ID)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "Spatial25 Cluster Proportions by arb_ID", x = "Spatial Cluster", y = "Proportion")

dev.off()